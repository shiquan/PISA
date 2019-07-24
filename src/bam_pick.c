// pick alignment reads by cell barcodes
#include "utils.h"
#include "number.h"
#include "thread_pool.h"
#include "barcode_list.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/khash_str2int.h"
#include "htslib/kseq.h"

static int usage()
{
    fprintf(stderr, "PickBam -list barcode.txt -tag CB -o out.bam in.bam\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *barcode_fname;
    const char *tag;

} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .barcode_fname = NULL,
    .tag = NULL,
};
static int parse_args(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-list") == 0) var = &args.barcode_fname;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
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

    // CHECK_EMPTY(args.barcode_fname, "-list Cell barcode list must be set.");
    CHECK_EMPTY(args.input_fname, "Input BAM file must be set.");
    CHECK_EMPTY(args.output_fname, "Output BAM file must be set.");
    
    return 0;
}

int bam_pick(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return usage();

    htsFile *fp = hts_open(args.input_fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");
    
    htsFile *out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(out, "%s : %s.", args.output_fname, strerror(errno));
    if (sam_hdr_write(out, hdr)== -1) error("Failed to write SAM header.");

    struct barcode_list *barcode = NULL;
    if (args.barcode_fname) {
        barcode = barcode_init();
        if (barcode_read(barcode, args.barcode_fname)) error("Empty barcode");
    }
    
    bam1_t *b;
    bam1_core_t *c;
    int ret;
    b = bam_init1();
    c = &b->core;
    int id;
    while ((ret = sam_read1(fp, hdr, b)) >= 0) {
        uint8_t *tag = bam_aux_get(b, args.tag);
        // if (!tag) error("Tag %s not found at line %s:%d\n", args.tag, hdr->target_name[c->tid], c->pos+1);
        if (barcode) {
            id = barcode_select(barcode, (char*)(tag+1));
            if (id == -1) continue;
            if (sam_write1(out, hdr, b) == -1) error("Failed to write SAM.");
        }
        else if (tag) {
            if (sam_write1(out, hdr, b) == -1) error("Failed to write SAM.");
        }
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    sam_close(out);
    if (barcode)
        barcode_destory(barcode);
    
    if (ret != -1) warnings("Truncated file?");
    return 0;
}
