#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "number.h"
#include <sys/stat.h>
#include "dict.h"
#include "pisa_version.h"

struct fp_out {
    BGZF *barcode_fp;
    BGZF *feature_fp;
    BGZF *mex_fp;
};

struct files {
    char *path;
    char *bc;
};

static int bins[4] = {50000, 100000, 500000, 1000000};

struct args {
    const char *input_fname;
    const char *bin_str;
    const char *outdir;
    const char *sample_list;
    const char *umi_tag;
    const char *bc_tag;
    const char *bc_list;
    int n_file;
    struct files *files;
    int mapq_thres;
    int n_bin;
    int *bins;
    struct fp_out *fps;
    int *n_records;
    int *n_features;
    int strand;
    int n_thread;
    struct dict *barcodes;
    int fix_barcodes;
} args = {
    .input_fname = NULL,
    .bin_str = NULL,
    .outdir = NULL,
    .sample_list = NULL,
    .umi_tag = NULL,
    .bc_tag = NULL,
    .bc_list = NULL,
    .n_file = 0,
    .mapq_thres = 10,
    .n_bin = 4,
    .bins = bins,
    .fps = NULL,
    .n_records = NULL,
    .n_features = NULL,
    .strand = 0,
    .n_thread = 10,
    //.hdr = NULL,
    .barcodes = NULL,
    .fix_barcodes = 0,
    // .in = NULL,
    //.idx = NULL,
};

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    const char *threads = NULL;
    const char *bin_str = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        
        if (strcmp(a, "-list") == 0) var = &args.bc_list;
        else if (strcmp(a, "-sample-list") == 0) var = &args.sample_list;
        else if (strcmp(a, "-cb") == 0) var = &args.bc_tag;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-bins") == 0) var = &bin_str;
        else if (strcmp(a, "-outdir") == 0) var = &args.outdir;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-@") == 0) var = &threads;
        else if (strcmp(a, "-strand") == 0) {
            args.strand = 1;
            continue;
        }

        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        
        error("Unknown argument, %s", a);
    }

    if (args.input_fname == NULL && args.sample_list == NULL)
        error("No input bam.");
    
    if (threads) args.n_thread = str2int((char*)threads);
    
    if (mapq) {
        args.mapq_thres = str2int(mapq); 
    }

    if (bin_str) {
        kstring_t tmp = {0,0,0};
        kputs(bin_str, &tmp);
        int n;
        int *s = ksplit(&tmp, ',', &n);
        if (n > 0 && n <= 4) {
            int k;
            for (k = 0; k < n; ++k) {
                char *ss = tmp.s + s[k];
                args.bins[k] = genome2int(ss);
            }
            args.n_bin = n;
        } else if (n > 4) {
            int *bins = malloc(n*sizeof(int));
            int k;
            for (k = 0; k < n; ++k) {
                char *ss = tmp.s + s[k];
                bins[k] = genome2int(ss);
            }
            args.bins = bins;
            args.n_bin = n;
        }
        free(s);
        free(tmp.s);
    }
    
    args.barcodes = dict_init();
    dict_set_value(args.barcodes);

    if (args.bc_tag) {
        if (args.bc_list) {
            int val = 0;
            dict_read2(args.barcodes, args.bc_list, &val);
            if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");
            args.fix_barcodes = 1;
        }
    } else {
        if (args.sample_list == NULL) {
            // no cell barcode, one sample
            dict_push(args.barcodes,"sample"); // make a pesudo sample name
        }
    }

    if (args.sample_list) {
        
    }


    if (args.input_fname) {
    }

    if (args.outdir) {
        struct stat sb;
        if (stat(args.outdir, &sb) != 0) error("Directory %s is not exist.", args.outdir);
        if (S_ISDIR(sb.st_mode) == 0) error("%s does not look like a directory.", args.outdir);
    }

    args.n_records = malloc(args.n_bin *sizeof(int));
    args.n_features = malloc(args.n_bin *sizeof(int));
    int k;
    for (k = 0; k < args.n_bin; ++k) {
        args.n_records[k] = 0;
        args.n_features[k] = 0;
    }
    
    args.fps = malloc(args.n_bin* sizeof(struct fp_out));

    kstring_t bc_str = {0,0,0};
    kstring_t ft_str = {0,0,0};
    kstring_t mx_str = {0,0,0};
    
    for (k = 0; k < args.n_bin; ++k) {
        bc_str.l = ft_str.l = mx_str.l = 0;
        if (args.outdir) {
            kputs(args.outdir, &bc_str);
            kputs(args.outdir, &ft_str);
            kputs(args.outdir, &mx_str);
            if (bc_str.s[bc_str.l-1] != '/') {
                kputc('/', &bc_str);
                kputc('/', &ft_str);
                kputc('/', &mx_str);
            }
        }
        kputw(args.bins[k], &bc_str);
        kputw(args.bins[k], &ft_str);
        kputw(args.bins[k], &mx_str);

        struct stat sb;
        if (stat(bc_str.s, &sb) != 0) {
            if (mkdir(bc_str.s, 0755)) error("%s : %s.", bc_str.s, strerror(errno));
        } else {
            if (S_ISDIR(sb.st_mode) == 0) error("%s does not look like a directory.", bc_str.s);
        }
        
        kputs("/barcodes.tsv.gz", &bc_str);
        kputs("/features.tsv.gz", &ft_str);
        kputs("/matrix.mtx.gz", &mx_str);
        
        struct fp_out *fps = &args.fps[k];
        fps->barcode_fp = bgzf_open(bc_str.s, "w");
        fps->feature_fp = bgzf_open(ft_str.s, "w");
        fps->mex_fp = bgzf_open(mx_str.s, "w");
    }

    free(bc_str.s);
    free(ft_str.s);
    free(mx_str.s);
    
    return 0;
}

static void memory_release()
{
    if (args.n_bin>4) free(args.bins);
    // close fps
    int i;
    for (i = 0; i < args.n_bin; ++i) {
        struct fp_out *fps = &args.fps[i];
        bgzf_close(fps->barcode_fp);
        bgzf_close(fps->feature_fp);
        bgzf_close(fps->mex_fp);
    }
    free(args.fps);
    free(args.n_records);
    free(args.n_features);

    //bam_hdr_destroy(args.hdr);
    // hts_close(args.in);

    for (i = 0; i < dict_size(args.barcodes); ++i) {
        struct cnt0 *smp = dict_query_value(args.barcodes, i);
        if (smp == NULL) continue;
        free(smp);
    }
    dict_destroy(args.barcodes);
}

struct cnt0 {
    int id; // cell id or sample id
    int cnt;
    struct dict *umi; // only active when umi set
};

struct cnt {
    int n, m;
    struct cnt0 *cnt;
};

struct chrom_counts {
    int need_update;
    int size; // chrom size, used for calculating bin number
    struct cnt **cnt; // windows, bins
};

static int cmpfunc(const void *_a, const void *_b)
{
    const struct cnt0 *a = (const struct cnt0*) _a;
    const struct cnt0 *b = (const struct cnt0*) _b;
    return (a->id > b->id) - (a->id < b->id);
}
struct chrom_counts *new_cc(int size)//, int id)
{
    struct chrom_counts *cc = malloc(sizeof(*cc));
    memset(cc, 0, sizeof(*cc));
    cc->size = size;
    cc->need_update = 0;
    cc->cnt = malloc(sizeof(struct cnt*)*args.n_bin);
    
    int i;
    for (i = 0; i < args.n_bin; ++i) {
        int bin = args.bins[i];
        int size0 = size/bin +1;
        args.n_features[i] += size0;
        cc->cnt[i] = malloc(sizeof(struct cnt)*size0);
        
        int j;
        for (j = 0; j < size0; ++j) {
            struct cnt *cnt = &cc->cnt[i][j];
            cnt->n = 0;
            cnt->m = 0;
            cnt->cnt = NULL;
        }
    }
    return cc;
}

void update_cc(struct chrom_counts *cc)
{
    if (cc->need_update) {
        int i;
        for (i = 0; i < args.n_bin; ++i) {
            int bin = args.bins[i];
            int size0 = cc->size/bin +1;
            int j;
            for (j = 0; j < size0; ++j) {
                struct cnt *cnt = &cc->cnt[i][j];
                int k;
                for (k = 0; k < cnt->n; ++k) {
                    struct cnt0 *cnt0 = &cnt->cnt[k];
                    if (cnt0->umi) {
                        cnt0->cnt = dict_size(cnt0->umi);
                        dict_destroy(cnt0->umi);
                        cnt0->umi = NULL;
                    }
                }
            }
        }
    }
}
void destroy_cc(struct chrom_counts *cc)
{
    int i;
    for (i = 0; i < args.n_bin; ++i) {
        int bin = args.bins[i];
        int size0 = cc->size/bin +1;
        int j;
        for (j = 0; j < size0; ++j) {
            struct cnt *cnt = &cc->cnt[i][j];
            if (cnt->m > 0) {
                free(cnt->cnt);
            }
        }
        free(cc->cnt[i]);
    }
    free(cc->cnt);
    free(cc);
    cc = NULL;
}
void write_out(bam_hdr_t *hdr, struct chrom_counts **ccs, int nref)
{
    kstring_t tmp = {0,0,0};
    if (args.outdir) {
        kputs(args.outdir, &tmp);
        if (tmp.s[tmp.l-1] != '/') kputc('/', &tmp);
    }
    kputs("meta.tsv.gz", &tmp);
    BGZF *fp = bgzf_open(tmp.s, "w");

    kstring_t bc_str = {0,0,0};
    // write count per cell/sample
    int i;
    tmp.l = 0;
    for (i = 0; i < dict_size(args.barcodes); ++i) {
        char *name = dict_name(args.barcodes, i);
        kputs(name, &bc_str);
        kputc('\n', &bc_str);

        struct cnt0 *smp = dict_query_value(args.barcodes, i);
        if (smp == NULL) continue;
        if (smp->umi && smp->cnt == 0) {
            smp->cnt = dict_size(smp->umi);
            dict_destroy(smp->umi);
        }
        kputs(name, &tmp);
        kputc('\t', &tmp);
        kputw(smp->cnt, &tmp);
        kputc('\n', &tmp);
    }
    kputs("", &tmp);
    kputs("", &bc_str);
    
    int l = bgzf_write(fp, tmp.s, tmp.l);
    if (l != tmp.l) error("Failed to write.");
    bgzf_close(fp);
    free(tmp.s);
    
    kstring_t ft_str = {0,0,0};
    kstring_t mx_str = {0,0,0};
    int k;
    for (k = 0; k < args.n_bin; ++k) {        
        ft_str.l = mx_str.l = 0;
        struct fp_out *fps = &args.fps[k];
        int bin = args.bins[k];
        // header

        kputs("%%MatrixMarket matrix coordinate integer general\n", &mx_str);
        kputs("% Generated by PISA ", &mx_str);
        kputs(PISA_VERSION, &mx_str);
        kputc('\n', &mx_str);
        ksprintf(&mx_str, "%d\t%d\t%d\n", args.n_features[k], dict_size(args.barcodes), args.n_records[k]);
        
        for (i = 0; i < nref; ++i) {
            struct chrom_counts *cc = ccs[i];
            if (cc == NULL) continue;
            int size = cc->size;
            int n = size/bin + 1;
            // struct cnt *cnt = cc->cnt[k];
            int j;
            for (j = 0; j < n; ++j) {
                int end = bin*(j+1);
                if (end > size) end = size;
                    
                if (args.strand) { //todo
                    ksprintf(&ft_str, "%s:%d-%d/+\n", hdr->target_name[i], bin*j+1, end);
                } else {
                    ksprintf(&ft_str, "%s:%d-%d\n", hdr->target_name[i], bin*j+1, end);
                }
                struct cnt *cnt0 = &cc->cnt[k][j];
                int j0;
                for (j0 = 0; j0 < cnt0->n; ++j0) {
                    struct cnt0 *cnt1 = &cnt0->cnt[j0];
                    ksprintf(&mx_str, "%d\t%d\t%d\n", j+1, cnt1->id+1, cnt1->cnt);
                }
            }            
        }

        l = bgzf_write(fps->feature_fp, ft_str.s, ft_str.l);
        if (l != ft_str.l) error("Failed to write file.");
        
        l = bgzf_write(fps->mex_fp, mx_str.s, mx_str.l);
        if (l != mx_str.l) error("Failed to write file.");

        l = bgzf_write(fps->barcode_fp, bc_str.s, bc_str.l);
        if (l != bc_str.l) error("Failed to write file.");                
    }

    if (bc_str.m) free(bc_str.s);
    if (ft_str.m) free(ft_str.s);
    if (mx_str.m) free(mx_str.s);

    for (i = 0; i > nref; ++i) {
        struct chrom_counts *cc = ccs[i];
        if (cc == NULL) continue;
        destroy_cc(cc);
    }
    free(ccs);
}

extern int bin_usage();
int bin_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bin_usage();

    htsFile *in = hts_open(args.input_fname, "r");
    if (in == NULL) error("%s : %s.", args.input_fname, strerror(errno));

    bam_hdr_t *hdr = sam_hdr_read(in);
    //hts_idx_t *idx = sam_index_load(in, args.input_fname);
    //if (idx == NULL) error("Failed to load bam index of %s", args.input_fname);

    int nref = sam_hdr_nref(hdr);
    /* int ll; */
    /* for (ll = 0; ll < nref; ++ll) { */
    /*     debug_print("%s\t%u", hdr->target_name[ll], idx->lidx[ll].offset[0]); */
    /* } */

    int checked[nref];
    memset(checked, 0, nref*sizeof(int));
    
    int lst_chr = -1;
    int lst_pos = -1;
    struct chrom_counts **ccs = malloc(sizeof(struct chrom_counts*)*nref);
    memset(ccs, 0, sizeof(struct chrom_counts*)*nref);

    // int fi;
    // for (fi = 0; fi < args.n_file; ++fi) {
        
    // #pragma omp parallel num_threads(args.n_thread)
        
    struct chrom_counts *cc = NULL;

    bam1_t *b;
    b = bam_init1();
    
    for (;;) {
        if (sam_read1(in, hdr, b) < 0) break;
        bam1_core_t *c = &b->core;
        if (c->tid < 0 || c->tid > nref) continue;
        if (c->qual < args.mapq_thres) continue;
        if (c->flag & BAM_FSECONDARY) continue;
        if (c->flag & BAM_FUNMAP) continue;
        if (c->flag & BAM_FQCFAIL) continue;
        if (c->flag & BAM_FDUP) continue;

        if (lst_chr == -1) {
            lst_chr = c->tid;
            checked[lst_chr] = 1;
            cc = new_cc(hdr->target_len[c->tid]);
            ccs[lst_chr] = cc;
            if (args.umi_tag) cc->need_update = 1;
            debug_print("Working on %s.", hdr->target_name[lst_chr]);
        }
        
        if (lst_chr != c->tid) {
            if (checked[c->tid]) error("BAM file is not sorted!");
            update_cc(cc);
            lst_chr = c->tid;
            checked[lst_chr] = 1;
            cc = new_cc(hdr->target_len[c->tid]);
            ccs[lst_chr] = cc;
            if (args.umi_tag) cc->need_update = 1;
            debug_print("Working on %s.", hdr->target_name[lst_chr]);
        }

        int pos = c->pos + 1;
        // check sorted
        if (lst_pos > pos) error("BAM file is not sorted!");

        // todo: junction reads
        
        // check barcode
        int id = 0;
        char *umi = NULL;
        
        if (args.bc_tag) {
            uint8_t *data = bam_aux_get(b, args.bc_tag);
            if (data == NULL) continue;
            
            id = dict_query(args.barcodes, (char*)(data+1));
            if (id == -1) {
                if (args.fix_barcodes) continue;
                id = dict_push(args.barcodes, (char*)(data+1));
            }
        }
        
        // check umi
        if (args.umi_tag) {
            uint8_t *data = bam_aux_get(b, args.umi_tag);
            if (data == NULL) continue;
            umi = (char*)(data+1);
        }

        struct cnt0 *smp = dict_query_value(args.barcodes, id);
        if (smp == NULL) {
            smp = malloc(sizeof(struct cnt0));
            memset(smp, 0, sizeof(struct cnt0));
            dict_assign_value(args.barcodes, id, smp);
        }
        
        smp->id = id;
        if (umi) {
            if (smp->umi == NULL) smp->umi = dict_init();
            dict_push(smp->umi, umi);
        } else {
            smp->cnt++;
        }
        
        int j;
        for (j = 0; j < args.n_bin; ++j) {
            int bin = pos/args.bins[j];
            struct cnt *cnt = &cc->cnt[j][bin];
            if (cnt->n == cnt->m) {
                if (cnt->m == 0) {
                    cnt->m = 2;
                    cnt->cnt = malloc(sizeof(struct cnt0)*cnt->m);
                } else {
                    cnt->m = cnt->m *2;
                    cnt->cnt = (struct cnt0*)realloc(cnt->cnt, sizeof(struct cnt0)*cnt->m);
                }
            }

            int k;
            for (k = 0; k < cnt->n; ++k) {
                struct cnt0 *cnt0 = &cnt->cnt[k];
                if (cnt0->id == id) {
                    if (umi) {
                        dict_push(cnt0->umi, umi);
                    } else {
                        cnt0->cnt++;
                    }
                    break;
                } else if (id < cnt0->id) {
                    struct cnt0 *cnt0 = &cnt->cnt[cnt->n++];
                    // thread safe
                    args.n_records[j]++;
                    //if (cnt->n==1) args.n_features[j]++;
                    
                    cnt0->cnt = 0;
                    cnt0->umi = NULL;
                    if (umi) {
                        cnt0->umi = dict_init();
                        dict_push(cnt0->umi, umi);
                    } else {
                        cnt0->cnt++;
                    }
                    cnt0->id = id;
                    
                    qsort(cnt->cnt, cnt->n, sizeof(struct cnt0), cmpfunc);
                    break;
                }
            }

            if (k == cnt->n) {
                if (cnt->n == cnt->m) {
                    cnt->m = cnt->m == 0 ? 2 : cnt->m *2;
                    cnt->cnt = realloc(cnt->cnt, sizeof(struct cnt0)*cnt->m);
                }

                struct cnt0 *cnt0 = &cnt->cnt[cnt->n++];
                args.n_records[j]++;
                //if (cnt->n==1) args.n_features[j]++;
                
                cnt0->cnt = 0;
                cnt0->umi = NULL;
                if (umi) {
                    cnt0->umi = dict_init();
                    dict_push(cnt0->umi, umi);
                } else {
                    cnt0->cnt++;
                }
                cnt0->id = id;
            } 

            
        }
        // fusion
        
    }
    
    update_cc(cc);

    write_out(hdr, ccs, nref); // also free ccs inside
    
    bam_destroy1(b);

    hts_close(in);
    bam_hdr_destroy(hdr);

    memory_release();

    return 0;
}
