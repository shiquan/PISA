#include "utils.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "dict.h"
#include "number.h"

static int usage()
{
    fprintf(stderr, "* Convert sam record to fragment file.");
    fprintf(stderr, "bam2frag in.bam\n");
    fprintf(stderr, "    -o        Output file. This file will be bgzipped and indexed.\n");
    fprintf(stderr, "    -S        Treat all input as single end.\n");
    fprintf(stderr, "    -tag      Cell barcode tag.\n");
    fprintf(stderr, "    -list     Cell barcode white list.\n");
    fprintf(stderr, "    -isize    Cutoff insert size to define NDRs.\n");
    fprintf(stderr, "    -@        Thread to unpack and pack files.\n");
    return 1;
}

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *barcode_list;
    const char *tag;
    int force_SE;
    int isize;
    int file_th;
    int qual_thres;
} args = {
    .input_fname    = NULL,
    .output_fname   = NULL,
    .barcode_list   = NULL,
    .tag            = "CB",
    .force_SE       = 0,
    .isize          = 0,
    .file_th        = 4,
    .qual_thres     = 20,
};

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *file_th = NULL;
    const char *isize   = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-list") == 0) var = &args.barcode_list;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-isize") == 0) var = &isize;
        else if (strcmp(a, "-S") == 0) {
            args.force_SE = 1;
            continue;
        }
        if (var != 0) {
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        
        error("Unknown argument : %s", a);
    }

    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (args.tag == NULL) error("-tag must be set.");
    if (strlen(args.tag) != 2) error("Bad format of tag, %s", args.tag);
    if (file_th) args.file_th = str2int(file_th);
    if (isize) args.isize = str2int(isize);
    if (args.isize < 0) args.isize = 0;
    
    return 0;
}  

struct frag {
    int id;
    int start;
    int end;
    int CB;
};

struct frag_pool {
    int n, m;
    struct frag *bed;
    struct dict *CB_dict;
};

static void frag_pool_push(struct frag_pool *p, bam1_t *b, int CB)
{
    if (p->n == p->m) {
        p->m = p->m == 0 ? 1024 : p->m<<1;
        p->bed = realloc(p->bed, p->m*sizeof(struct frag));
    }
    struct frag *bed = &p->bed[p->n];
    bed->id = b->core.tid;
    bed->start = b->core.pos;
    bed->end = b->core.isize == 0 ? bed->start + b->core.l_qseq : bed->start + b->core.isize;
    bed->CB = CB;
    p->n++;
}

static int cmpfunc(const void *_a, const void *_b)
{
    struct frag *a = (struct frag*)_a;
    struct frag *b = (struct frag*)_b;
    return a->start - b->start == 0 ? a->end - b->end : a->start - b->start;
}
static int frag_pool_print(BGZF *fp, struct frag_pool *p, bam_hdr_t *hdr)
{
    if (p->n == 0) return 1;
    qsort(p->bed, p->n, sizeof(struct frag), cmpfunc);
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < p->n; ++i) {
        str.l = 0;
        struct frag *bed = &p->bed[i];
        kputs(hdr->target_name[bed->id], &str); kputc('\t', &str);        
        kputw(bed->start, &str); kputc('\t', &str);
        kputw(bed->end, &str); kputc('\t', &str);
        kputs(dict_name(p->CB_dict, bed->CB), &str);
        kputs("\t1\n", &str);
        if (bgzf_write(fp, str.s, str.l) < 0) error("Failed to write file.");
    }
    free(str.s);
    p->n = 0;
    return 0;
}

int bam2frag(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    htsFile *fp  = hts_open(args.input_fname, "r");
    if (fp == NULL)
        error("%s : %s.", args.input_fname, strerror(errno));
    
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");

    hts_set_threads(fp, args.file_th);

    BGZF *fp_out = bgzf_open(args.output_fname, "w");
    if (fp_out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    
    bgzf_mt(fp_out, args.file_th, 256);

    bam1_t *b = bam_init1();

    struct frag_pool *p = malloc(sizeof(struct frag_pool));
    memset(p, 0, sizeof(*p));
    p->CB_dict = dict_init();
    int white_list_flag = 0;
    if (args.barcode_list) {
        dict_read(p->CB_dict, args.barcode_list);
        if (dict_size(p->CB_dict) == 0) error("Empty barcode list.");
        white_list_flag = 1;
    }
    int last_id = -1;
    int ret;
    while ((ret = sam_read1(fp, hdr, b)) >=0) {
        if (args.qual_thres > 0 && b->core.qual < args.qual_thres) continue;
        if (b->core.tid < 0) continue;
        if (args.force_SE == 0) {
            if (args.isize && b->core.isize > args.isize) continue;
            if (b->core.isize < 0) continue;
        }
        if (last_id == -1) last_id = b->core.tid;
        // output buffered Records
        if (last_id != b->core.tid) {
            frag_pool_print(fp_out, p, hdr);
            last_id = b->core.tid;
        }
        
        uint8_t *tag = bam_aux_get(b, args.tag);
        if (!tag) continue;
        int cb = dict_query(p->CB_dict, (char*)(tag+1));
        if (cb == -1) {
            if (white_list_flag) continue;
            cb = dict_push(p->CB_dict, (char*)(tag+1));
        }
        
        frag_pool_push(p, b, cb);
    }    
    frag_pool_print(fp_out,p,hdr);
    bgzf_close(fp_out);        
    bam_destroy1(b);
    dict_destroy(p->CB_dict);
    free(p->bed);
    free(p);
    bam_hdr_destroy(hdr);
    hts_close(fp);

    const tbx_conf_t tbx_conf_bed = { TBX_UCSC, 1, 2, 3, '#', 0 };
    
    if (tbx_index_build(args.output_fname, 0, &tbx_conf_bed))
        warnings("Failed to build index file of %s.", args.output_fname);
    
    return 0;
}
