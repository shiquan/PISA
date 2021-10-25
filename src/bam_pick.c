// pick alignment reads by cell barcodes
#include "utils.h"
#include "number.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/khash_str2int.h"
#include "htslib/kseq.h"
#include "dict.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8193);

static struct args {
    const char *input_fname;
    const char *output_fname;

    char **tags;
    int n_tag;
    struct dict *barcodes;
    int file_th;

    htsFile *fp_in;
    htsFile *fp_out;
    bam_hdr_t *hdr;

    int mapq_thres;

    int check_wl_column; 
} args = {
    .input_fname = NULL,
    .output_fname = NULL,

    .tags = NULL,
    .n_tag = 0,
    .file_th = 4,

    .barcodes = NULL,
    .fp_in = NULL,
    .fp_out = NULL,
    .hdr = NULL,
    .mapq_thres = 0,
    .check_wl_column = 0,
};

extern int pick_usage();

static int parse_args(int argc, char **argv)
{
    const char *file_th = NULL;
    const char *tag_str = NULL;
    const char *extract_fname = NULL;
    const char *mapq = NULL;
    int i;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-list") == 0) var = &extract_fname;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-tags") == 0) var = &tag_str;
        else if (strcmp(a, "-tag") == 0) error("Did you mean -tags ?");
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument: %s", a);
    }

    CHECK_EMPTY(args.input_fname, "Input BAM file must be set.");
    CHECK_EMPTY(args.output_fname, "Output BAM file must be set.");

    if (tag_str == NULL) error("-tag is required.");

    if (file_th) args.file_th = str2int((char*)file_th);
    if (mapq) args.mapq_thres = str2int(mapq);
    
    args.fp_in = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp_in, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp_in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    args.fp_out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.fp_out, "%s : %s.", args.output_fname, strerror(errno));

    hts_set_threads(args.fp_in, args.file_th);

    args.hdr = sam_hdr_read(args.fp_in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");

    if (sam_hdr_write(args.fp_out, args.hdr)) error("Failed to write SAM header.");

    if (extract_fname) {
        gzFile fp = gzopen(extract_fname, "r");
        if (fp == NULL) error("%s : %s.", extract_fname, strerror(errno));
        
        kstream_t *ks = ks_init(fp);
        
        kstring_t str = {0,0,0};
        kputs(tag_str, &str);
        
        int n;
        int *s = ksplit(&str, ',', &n);
        
        args.n_tag = n;
        args.tags = malloc(n*sizeof(void**));
        
        for (i = 0; i < n; ++i) {
            args.tags[i] = strdup(str.s+s[i]);
        }
        
        free(s);
        
        args.barcodes = dict_init();
        
        int ret;
        kstring_t tmp = {0,0,0};
        while (ks_getuntil(ks, 2, &str, &ret) >= 0) {
            int *s = ksplit(&str, '\t', &n);
            if (n < args.n_tag) {
                args.check_wl_column = n;
            }
            else {
                args.check_wl_column = args.n_tag;
            }
                //error("Malformed line?");
            
            tmp.l = 0;
            
            int i;
            for (i = 0; i < n; ++i) {
                kputs(str.s+s[i], &tmp);
            }
            dict_push(args.barcodes, tmp.s);
        }
        if (dict_size(args.barcodes) == 0) error("Empty list.");
        
        free(tmp.s);
        free(str.s);
        gzclose(fp);
        ks_destroy(ks);
    }

    return 0;
}
static void memory_release()
{
    if (args.barcodes)
        dict_destroy(args.barcodes);
    bam_hdr_destroy(args.hdr);
    sam_close(args.fp_in);
    sam_close(args.fp_out);
}

int bam_pick(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return pick_usage();

    kstring_t str = {0,0,0};
    
    bam1_t *b;
    b = bam_init1();
    int ret;
    int skip_flag;
    for (;;) {
        // init statue
        str.l = 0;
        skip_flag = 0;
        
        ret = sam_read1(args.fp_in, args.hdr, b);
        if (ret < 0) break;
        
        if (b->core.flag & BAM_FQCFAIL) continue;
        if (b->core.flag & BAM_FSECONDARY) continue;
        if (b->core.qual < args.mapq_thres) continue;
        int i;
        //for (i = 0; i < args.n_tag; ++i) {
        for (i = 0; i < args.check_wl_column; ++i) {
            uint8_t *tag = bam_aux_get(b, args.tags[i]);
            if (!tag) {
                str.l = 0;
                skip_flag = 1;
                break;
            }
            if (*tag == 'A') kputc(tag[1], &str);
            else kputs((char*)(tag+1), &str);
        }

        if (skip_flag == 1) continue;
        // check tag exists
        for ( ; i < args.n_tag; ++i) {
            uint8_t *tag = bam_aux_get(b, args.tags[i]);
            if (!tag) {
                str.l = 0;
                skip_flag = 1;
                break;
            }
        }

        if (skip_flag == 1) continue;
        if(str.l) {
            // query barcodes
            int query = dict_query(args.barcodes, str.s);
            if (query < 0) continue; // jump to next record
        }
        // output
        if (sam_write1(args.fp_out, args.hdr, b) == -1) error("Failed to write SAM.");      
    }
    bam_destroy1(b);
    free(str.s);
    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());

    return 0;
}
