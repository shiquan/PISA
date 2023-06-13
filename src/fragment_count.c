#include "utils.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "dict.h"
#include "bed.h"
#include "region_index.h"
#include "number.h"
#include "pisa_version.h"
#include <omp.h>

static struct args {
    const char *input_fname;
    const char *bed_fname;
    const char *barcode_list;
    const char *outdir;
    const char *prefix;
    struct bed_spec *B;
    int n_thread;
} args = {
    .input_fname  = NULL,
    .bed_fname    = NULL,
    .barcode_list = NULL,
    .outdir       = NULL,
    .prefix       = NULL,
    .B            = NULL,
    .n_thread     = 4,
};

static int wl = 0;

static int parse_args(int argc, char **argv)
{
    int i;
    const char *threads = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-list") == 0) var = &args.barcode_list;
        else if (strcmp(a, "-outdir") == 0) var = &args.outdir;
        else if (strcmp(a, "-prefix") == 0) var = &args.prefix;
        else if (strcmp(a, "-t") == 0) var = &threads;
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == 0) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument, %s.", a);
    }

    if (args.input_fname == NULL) error("No input fragment specified.");

    if (threads) args.n_thread = str2int(threads);
    // check index
    tbx_t *tbx = tbx_index_load(args.input_fname);
    if (tbx == NULL) {
        warnings("Fail to load index file, try to index it ..");
        if (tbx_index_build(args.input_fname, 0, &tbx_conf_bed))
            error("Failed to index.");
    } else {
        tbx_destroy(tbx);
    }
    if (args.bed_fname == NULL) error("No -bed specified.");
    args.B = bed_read(args.bed_fname);
        
    if (args.barcode_list) wl = 1;
    
    return 0;
}

extern int frag_count_usage();

struct ret {
    uint32_t counts;
    struct dict *bc;
};

struct ret *query_count_write(struct bed_spec *B, int ci, const char *fn)
{
    kstring_t str = {0,0,0};
    kstring_t r = {0,0,0};
    BGZF *fp = bgzf_open(args.input_fname, "r");
    if (fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
            
    tbx_t *tbx = tbx_index_load(args.input_fname);
    if (tbx == NULL) error("Failed to load index file.");
    
    struct dict *cc = dict_init();
    dict_set_value(cc);
    
    BGZF *out = bgzf_open(fn, "w");
    if (out == NULL) error("%s : %s.", fn, strerror(errno));

    struct ret *ret = malloc(sizeof(struct ret));
    ret->counts = 0;
    ret->bc = dict_init();

    if (wl == 1) dict_read(cc, args.barcode_list, 0);
    
    int j;
    j = B->ctg[ci].idx ;
    int i = 0;

    while (1) {
        
        if (i == B->ctg[ci].offset) break;
        const struct bed *b = &B->bed[j++];
        i++;
        char *name = bed_seqname(B, b->seqname);
        int id = tbx_name2id(tbx, name);
        hts_itr_t *itr = tbx_itr_queryi(tbx, id, b->start, b->end);
            
        while(1) {
            int ret = tbx_bgzf_itr_next(fp, tbx, itr, &r);
            if (ret < 0) break;
            int n;
            int *s = ksplit(&r, '\t', &n);
            assert(n >= 4);
            int start = str2int(r.s + s[1]);
            int end = str2int(r.s + s[2]);
            
            if (end <= b->start) continue;
            if (start >= b->end) continue;
            
            char *cell = r.s + s[3];
            
            int idx = dict_query(cc, cell);
            free(s);
            
            if (idx < 0) {
                if (wl == 0) idx = dict_push(cc, cell);
                else continue;
            }
            
            int *v = dict_query_value(cc, idx);
            if (v == NULL) {
                v = malloc(sizeof(int));
                *v = 1;
                dict_assign_value(cc, idx, v);
            } else {
                *v = *v + 1;
            }
        }
        
        // print out, init value to 0
        int k = 0;

        while (1) {
            if (k == dict_size(cc)) break;
            int *v = dict_query_value(cc, k);
            if (v != NULL && *v != 0) {
                ksprintf(&str,"%d\t%s\t%d\n", j , dict_name(cc,k), *v); // j already increased
                ret->counts ++;
                *v = 0;
            }
            ++k; // increase iterator
        }
        int size = bgzf_write(out, str.s, str.l);
        if (size != str.l) warnings("Write size is wrong!");
        
        str.l = 0;
        hts_itr_destroy(itr);
    }
            
    if (r.m) free(r.s);
    if (str.m) free(str.s);
    bgzf_close(fp);
    bgzf_close(out);
    tbx_destroy(tbx);

    for (i = 0; i < dict_size(cc); ++i) {
        int *v = dict_query_value(cc, i);
        char *name = dict_name(cc,i);
        dict_push(ret->bc, name);
        if (v) free(v);
    }
    dict_destroy(cc);
    return ret;
}


struct ret *create_temp(struct bed_spec *B, const char **tmpfiles)
{
    int ci;
    int n = dict_size(B->seqname);

    struct ret *ret = malloc(sizeof(struct ret));
    ret->counts = 0;
    ret->bc = dict_init();
    
    omp_lock_t writelock;
    omp_init_lock(&writelock);

#pragma omp parallel for num_threads(args.n_thread) schedule(dynamic)
    for (ci = 0; ci < n; ++ci) {
        struct ret *ret0 = query_count_write(B, ci, tmpfiles[ci]);
        omp_set_lock(&writelock);
        ret->counts += ret0->counts;
        int j = 0;
        while(1) {
            if (j == dict_size(ret0->bc)) break;
            dict_push(ret->bc, dict_name(ret0->bc, j));
            j++;
        }
        dict_destroy(ret0->bc);
        free(ret0);
        LOG_print("Processed %s.", dict_name(B->seqname,ci));
        omp_unset_lock(&writelock);
    }

    omp_destroy_lock(&writelock);
    return ret;
}

int fragment_count(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return frag_count_usage();

    struct bed_spec *B = args.B;
    if (B->n ==0) error("Empty bed file.");

    kstring_t barcode_str = {0,0,0};
    kstring_t bed_str = {0,0,0};
    kstring_t mex_str = {0,0,0};
    kstring_t temp_str = {0,0,0};
    if (args.outdir) {
        kputs(args.outdir, &barcode_str);
        kputs(args.outdir, &bed_str);
        kputs(args.outdir, &mex_str);
        kputs(args.outdir, &temp_str);
        if (args.outdir[strlen(args.outdir)-1] != '/') {
            kputc('/', &barcode_str);
            kputc('/', &bed_str);
            kputc('/', &mex_str);
            kputc('/', &temp_str);
        }

        if (args.prefix) {
            kputs(args.prefix, &barcode_str);
            kputs(args.prefix, &bed_str);
            kputs(args.prefix, &mex_str);            
        }
    }

    kputs("barcodes.tsv.gz", &barcode_str);
    kputs("peaks.bed.gz", &bed_str);
    kputs("matrix.mtx.gz", &mex_str);
    kputs("__temp_", &temp_str);
    
    BGZF *fout_mex = bgzf_open(mex_str.s, "w");
    BGZF *fout_bc = bgzf_open(barcode_str.s, "w");
    BGZF *fout_bed = bgzf_open(bed_str.s, "w");

    if (fout_mex == NULL || fout_bc == NULL || fout_bed == NULL)
        error("Failed to create matrix files.");

    int n = dict_size(B->seqname);
    char **tmpfiles = malloc(sizeof(char**)*n);

    int i;
    for (i = 0; i < n; ++i) {
        kstring_t f = {0,0,0};
        kputs(temp_str.s, &f);
        kputw(i, &f);
        tmpfiles[i] = f.s;
    }

    struct ret *ret;
    ret = create_temp(B,(const char**)tmpfiles);

    int n_cell = dict_size(ret->bc);
    int n_feature = B->n;
    kstring_t str = {0,0,0};


    // write cells
    for (i = 0; i < n_cell; ++i) {
        kputs(dict_name(ret->bc,i), &str);
        kputc('\n', &str);
    }
    int size = bgzf_write(fout_bc, str.s, str.l);
    if (size != str.l) warnings("Size is wrong!");
    str.l = 0;

    bgzf_close(fout_bc);

    // write bed
    for (i = 0; i < B->n; ++i) {
        struct bed *b = &B->bed[i];
        char *name = bed_seqname(B, b->seqname);
        ksprintf(&str, "%s\t%d\t%d\n", name, b->start, b->end);        
    }
    size = bgzf_write(fout_bed, str.s, str.l);
    if (size != str.l) warnings("Size is wrong!");
    str.l = 0;
    bgzf_close(fout_bed);
    bgzf_mt(fout_mex, args.n_thread, 256);
    
    kputs("%%MatrixMarket matrix coordinate integer general\n", &str);
    kputs("% Generated by PISA ", &str);
    kputs(PISA_VERSION, &str);
    kputc('\n', &str);
    ksprintf(&str, "%d\t%d\t%u\n", n_feature, n_cell, ret->counts);

    size = bgzf_write(fout_mex, str.s, str.l);
    if (size != str.l) warnings("Size is wrong!");
    str.l = 0;

    LOG_print("Merge temp files ..");
    kstring_t temp = {0,0,0};
    for (i = 0; i < n; ++i) {
        int ncol;
        BGZF *fp = bgzf_open(tmpfiles[i], "r");
        while(1) {
            int _r = bgzf_getline(fp, '\n', &temp);
            if (_r <= 0) break;
            int *s = ksplit(&temp, '\t', &ncol);
            assert(ncol == 3);
            kputs(temp.s, &str);
            kputc('\t', &str);
            int id = dict_query(ret->bc,temp.s + s[1])+1;
            kputw(id, &str);
            kputc('\t', &str);
            kputs(temp.s+s[2], &str);
            kputc('\n', &str);
            
            size = bgzf_write(fout_mex, str.s, str.l);
            if (size != str.l) warnings("Size is wrong!");
            free(s);
            str.l = 0;
        }

        bgzf_close(fp);
        remove(tmpfiles[i]);
        free(tmpfiles[i]);
    }
    free(temp.s);
    free(str.s);
    free(tmpfiles);
    
    bgzf_close(fout_mex);

    dict_destroy(ret->bc);
    free(ret);
    bed_spec_destroy(B);

    free(barcode_str.s);
    free(bed_str.s);
    free(mex_str.s);
    free(temp_str.s);

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}

