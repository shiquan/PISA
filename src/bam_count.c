// count reads/fragments matrix for single-cell datasets
#include "utils.h"
#include "number.h"
#include "dict.h"
#include "dna_pool.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <sys/stat.h>
#include "pisa_version.h" // mex output

static struct args {
    const char *input_fname;
    const char *whitelist_fname;
    const char *output_fname;
    const char *outdir; // v0.4, support Market Exchange Format (MEX) for sparse matrices
        
    const char *tag; // cell barcode tag
    const char *anno_tag; // feature tag
    const char *umi_tag;

    struct dict *features;
    struct dict *barcodes;
    
    int mapq_thres;
    int use_dup;
    int enable_corr_umi;
    int n_thread;
    int one_hit;
    
    htsFile *fp_in;
    bam_hdr_t *hdr;

    uint64_t n_record;
} args = {
    .input_fname     = NULL,
    .whitelist_fname = NULL,
    .output_fname    = NULL,
    .outdir          = NULL,
    .tag             = NULL,
    .anno_tag        = NULL,
    .umi_tag         = NULL,

    .barcodes        = NULL,
    .features        = NULL,
    
    .mapq_thres      = 20,
    .use_dup         = 0,
    .enable_corr_umi = 0,
    .one_hit         = 0,
    .n_thread        = 5,
    .fp_in           = NULL,
    .hdr             = NULL,
    .n_record        = 0
};

union counts {
    uint32_t count;
    struct PISA_dna_pool *p;
};

static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.fp_in);
    
    int i;
    int n_feature;
    n_feature = dict_size(args.features);
    for (i = 0; i < n_feature; ++i) {
        struct PISA_dna_pool *v = dict_query_value(args.features, i);
        PISA_idx_destroy(v);
    }
    dict_destroy(args.features);
    dict_destroy(args.barcodes);
}
    
extern int bam_count_usage();

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    const char *n_thread = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-anno-tag") == 0) var = &args.anno_tag;
        else if (strcmp(a, "-list") == 0) var = &args.whitelist_fname;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-outdir") == 0) var = &args.outdir;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-@") == 0) var = &n_thread;
        else if (strcmp(a, "-dup") == 0) {
            args.use_dup = 1;
            continue;
        }
        else if (strcmp(a, "-one-hit") == 0) {
            args.one_hit = 1;
            continue;
        }
        
        else if (strcmp(a, "-corr") == 0) {
            //args.enable_corr_umi = 1;
            warnings("Option -corr has been removed since v0.8, to correct UMIs please use `PISA corr` instead.");
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
    if (args.input_fname == 0) error("No input bam.");
    if (args.output_fname) {
        warnings("PISA now support MEX format. Old cell X gene expression format is very poor performance. Try -outdir instead of -o.");
    }
    
    if (args.tag == 0) error("No cell barcode specified.");
    if (args.anno_tag == 0) error("No anno tag specified.");

    if (n_thread) args.n_thread = str2int((char*)n_thread);

    if (args.outdir) {
         struct stat sb;
         if (stat(args.outdir, &sb) != 0) error("Directory %s is not exists.", args.outdir);
         if (S_ISDIR(sb.st_mode) == 0)  error("%s does not look like a directory.", args.outdir);
    }
    
    args.fp_in = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp_in, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp_in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    hts_set_threads(args.fp_in, args.n_thread);
    
    args.hdr = sam_hdr_read(args.fp_in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");

    if (mapq) {
        args.mapq_thres = str2int(mapq);        
    }
    args.features = dict_init();
    dict_set_value(args.features);
    
    args.barcodes = dict_init();

    if (args.whitelist_fname) {
        dict_read(args.barcodes, args.whitelist_fname);
        if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");
    }
    
    return 0;
}

int *str_split(kstring_t *str, int *_n)
{
    int m=1, n=0;
    int i;
    int *s = malloc(1*sizeof(int));
    s[n++] = 0;
    for (i = 0; i < str->l; ++i) {
        if (str->s[i] == ',' || str->s[i] == ';') {
            if (m == n) {
                m += 2;
                s = realloc(s, sizeof(int)*m);
            }
            s[n++] = i+1;
            str->s[i] = '\0';
        }
    }
    *_n = n;
    return s;
}
int count_matrix_core(bam1_t *b)
{
    uint8_t *tag = bam_aux_get(b, args.tag);
    if (!tag) return 1;
        
    uint8_t *anno_tag = bam_aux_get(b, args.anno_tag);
    if (!anno_tag) return 1;

    if (args.umi_tag) {
        uint8_t *umi_tag = bam_aux_get(b, args.umi_tag);
        if (!umi_tag) return 1;
    }

    int cell_id;
    if (args.whitelist_fname) {
        cell_id = dict_query(args.barcodes, (char*)(tag+1));
        if (cell_id == -1) return 1;
    }
    else {
        cell_id = dict_push(args.barcodes, (char*)(tag+1));
    }

    // for each feature
    kstring_t str = {0,0,0};
    kputs((char*)(anno_tag+1), &str);
    int n_gene;
    int *s = str_split(&str, &n_gene); // seperator ; or ,

    // Sometime two or more genes or functional regions can overlapped with each other, if default PISA counts the reads for both of these regions.
    // But if -one-hit set, these reads will be filtered.
    if (args.one_hit == 1 && n_gene >1) {
        free(str.s);
        free(s);
        return 1;
    }
    
    int i;
    for (i = 0; i < n_gene; ++i) {
        // Features (Gene or Region)
        char *val = str.s + s[i];
        
        int idx = dict_query(args.features, val);
        if (idx == -1) idx = dict_push(args.features, val);

        struct PISA_dna_pool *v = dict_query_value(args.features, idx);

        if (v == NULL) {
            v = PISA_dna_pool_init();
            dict_assign_value(args.features, idx, v);
        } 
        // not store cell barcode for each hash, use id number instead to reduce memory
        struct PISA_dna *c= PISA_idx_query(v, cell_id);
        if (c == NULL) {
            c = PISA_idx_push(v, cell_id);
        //if (c->data == NULL) {
            union counts *counts = malloc(sizeof(union counts));
            if (args.umi_tag) 
                counts->p = PISA_dna_pool_init();
            else
                counts->count = 0;
            c->data = counts;
        }

        if (args.umi_tag) {
            uint8_t *umi_tag = bam_aux_get(b, args.umi_tag);
            assert(umi_tag);
            char *val = (char*)(umi_tag+1);
            assert(c->data);
            union counts *count = c->data;

            PISA_dna_push(count->p, val);
        }
        else {
            union counts *count = c->data;
            count->count++;
        }
    }
    free(str.s);
    free(s);
    return 0;
}

static void update_counts()
{
    int n_feature = dict_size(args.features);
    int i;
    for (i = 0; i < n_feature; ++i) {
        struct PISA_dna_pool *v = dict_query_value(args.features, i);
        int j;
        int n_cell = v->l;
        for (j = 0; j < n_cell; ++j) {
            union counts *count = v->data[j].data;
            assert(count);
            if (args.umi_tag) {
                int size = count->p->l;
                PISA_dna_destroy(count->p);
                count->count = size;
            }
        }
        args.n_record += n_cell;
    }
}
static void write_outs()
{
    int n_barcode = dict_size(args.barcodes);
    int n_feature = dict_size(args.features);

    if (n_barcode == 0) error("No barcode found.");
    if (n_feature == 0) error("No feature found.");
    if (args.n_record == 0) {
        warnings("No anntated record found.");
        return;
    }
    
    if (args.outdir) {
        kstring_t barcode_str = {0,0,0};
        kstring_t feature_str = {0,0,0};
        kstring_t mex_str = {0,0,0};
        kputs(args.outdir, &barcode_str);
        kputs(args.outdir, &feature_str);
        kputs(args.outdir, &mex_str);

        if (args.outdir[strlen(args.outdir)-1] != '/') {
            kputc('/', &barcode_str);
            kputc('/', &feature_str);
            kputc('/', &mex_str);
        }

        kputs("barcodes.tsv.gz", &barcode_str);
        kputs("features.tsv.gz", &feature_str);
        kputs("matrix.mtx.gz", &mex_str);
        
        BGZF *barcode_fp = bgzf_open(barcode_str.s, "w");
        bgzf_mt(barcode_fp, args.n_thread, 256);
        CHECK_EMPTY(barcode_fp, "%s : %s.", barcode_str.s, strerror(errno));
        
        int i;

        kstring_t str = {0,0,0};
        
        for (i = 0; i < n_barcode; ++i) {
            kputs(dict_name(args.barcodes, i), &str);
            kputc('\n', &str);
        }
        int l = bgzf_write(barcode_fp, str.s, str.l);
        if (l != str.l) error("Failed to write.");
        bgzf_close(barcode_fp);

        str.l = 0;
        BGZF *feature_fp = bgzf_open(feature_str.s, "w");
        bgzf_mt(feature_fp, args.n_thread, 256);
        CHECK_EMPTY(feature_fp, "%s : %s.", feature_str.s, strerror(errno));
        for (i = 0; i < n_feature; ++i) {
            kputs(dict_name(args.features,i), &str);
            kputc('\n', &str);
        }
        l = bgzf_write(feature_fp, str.s, str.l);
        if (l != str.l) error("Failed to write.");
        
        bgzf_close(feature_fp);

        str.l = 0;

        BGZF *mex_fp = bgzf_open(mex_str.s, "w");
        CHECK_EMPTY(mex_fp, "%s : %s.", mex_str.s, strerror(errno));
        
        bgzf_mt(mex_fp, args.n_thread, 256);
        kputs("%%MatrixMarket matrix coordinate integer general\n", &str);
        kputs("% Generated by PISA ", &str);
        kputs(PISA_VERSION, &str);
        kputc('\n', &str);
        ksprintf(&str, "%d\t%d\t%llu\n", n_feature, n_barcode, args.n_record);

        for (i = 0; i < n_feature; ++i) {
            struct PISA_dna_pool *v = dict_query_value(args.features, i);
            int j;
            int n_cell = v->l;
            for (j = 0; j < n_cell; ++j) {
                union counts *count = v->data[j].data;
                ksprintf(&str, "%d\t%d\t%u\n", i+1, v->data[j].idx+1, count->count);
            }

            if (str.l > 100000000) {
                int l = bgzf_write(mex_fp, str.s, str.l);
                if (l != str.l) error("Failed to write file.");
                str.l = 0;
            }
        }

        if (str.l) {
            l = bgzf_write(mex_fp, str.s, str.l);
            if (l != str.l) error("Failed to wirte.");
        }

        free(str.s);
        free(mex_str.s);
        free(barcode_str.s);
        free(feature_str.s);
        bgzf_close(mex_fp);
    }
    
    // header
    if (args.output_fname) {
        int i;
        FILE *out = fopen(args.output_fname, "w");
        CHECK_EMPTY(out, "%s : %s.", args.output_fname, strerror(errno));
        fputs("ID", out);
        
        for (i = 0; i < n_barcode; ++i)
            fprintf(out, "\t%s", dict_name(args.barcodes, i));
        fprintf(out, "\n");
        uint32_t *temp = malloc(n_barcode*sizeof(int));        
        for (i = 0; i < n_feature; ++i) {
            struct PISA_dna_pool *v = dict_query_value(args.features, i);
            int j;
            int n_cell = v->l;
            memset(temp, 0, sizeof(int)*n_barcode);
            fputs(dict_name(args.features, i), out);
            for (j = 0; j < n_cell; ++j) {
                int idx = v->data[j].idx;
                temp[idx] = v->data[j].count;
            }

            for (j = 0; j < n_barcode; ++j)
                fprintf(out, "\t%u", temp[j]);
            fputc('\n', out);
        }
        fclose(out);
        free(temp);
    }

}

int count_matrix(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return bam_count_usage();
        
    bam1_t *b;

    int ret;
    b = bam_init1();
    
    for (;;) {
        ret = sam_read1(args.fp_in, args.hdr, b);
        if (ret < 0) break;
                
        bam1_core_t *c;
        c = &b->core;

        if (c->tid <= -1 || c->tid > args.hdr->n_targets || (c->flag & BAM_FUNMAP)) continue;
        if (c->qual < args.mapq_thres) continue;
        if (args.use_dup == 0 && c->flag & BAM_FDUP) continue;
        count_matrix_core(b);
    }
    
    bam_destroy1(b);
    
    if (ret != -1) warnings("Truncated file?");   

    update_counts();

    write_outs();
    
    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
