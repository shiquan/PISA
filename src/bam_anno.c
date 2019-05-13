// annotate Gene or Peak to SAM attribution
#include "utils.h"
#include "number.h"
#include "thread_pool.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/khash_str2int.h"
#include "htslib/kseq.h"
#include "bed_lite.h"
#include <zlib.h>

static int usage()
{
    fprintf(stderr, "Annotate gene or peak information to SAM attribution\n");
    fprintf(stderr, "AnnoBam: -bed peak.bed -tag attr -o out.bam in.bam\n");
    fprintf(stderr, "  -o           Output BAM.");
    fprintf(stderr, "For peaks:\n");
    fprintf(stderr, "  -bed         Peak region or functional region in BED format.\n");
    fprintf(stderr, "  -tag         Tag name.\n");
    /*
    fprintf(stderr, "For genes:\n");
    fprintf(stderr, "  -gff         Genome annotation database in GFF format.\n");
    fprintf(stderr, "  -gtf         Genome annotation database in GTF format.\n");
    fprintf(stderr, "Notices:\n");
    fprintf(stderr, "* Default tags for genes: FC for functional region, GE for gene name.\n");
    fprintf(stderr, "* -bed, -gff, -gtf are conflict with each other.\n");
    */
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bed_fname;
    const char *tag; // attribute in BAM
    const char *gff_fname;
    const char *gtf_fname;
    struct bedaux *bed;
    // struct ganno *ganno;
    htsFile *fp;
    htsFile *out;
    bam_hdr_t *hdr;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .bed_fname = NULL,
    .tag = NULL,
    .gff_fname = NULL,
    .gtf_fname = NULL,

    .bed = NULL,
    // .ganno = NULL,
    .hdr = NULL,
};
static int parse_args(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-gff") == 0) var = &args.gff_fname;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
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
    
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    
    if (args.bed_fname == NULL && args.gff_fname == NULL && args.gtf_fname == NULL) error("One of -bed, -gtf and -gff must be set.");
    
    if (args.bed_fname) {        
        if (args.gff_fname || args.gtf_fname) error("-gff and -gtf is conflict with -bed.");
        CHECK_EMPTY(args.tag, "-tag must be specified if set -bed.");
        args.bed = bed_read(args.bed_fname);
        CHECK_EMPTY(args.bed, "Bed is empty.");
    }

    /*
    if (args.gtf_fname) {
        if (args.gff_fname) error("-gff is conflict with -gtf.");
         args.ganno = ganno_open_gtf(args.gtf_fname);
    }
    else
        args.ganno = ganno_open_gff(args.gff_fname);

    CHECK_EMPTY(args.ganno, "Failed to load annotation file.");
    */
    args.fp  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp);

    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    args.hdr = sam_hdr_read(args.fp);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
    
    return 0;
}
void bam_anno_attr_bed()
{
    bam1_t *b;
    bam1_core_t *c;
    int ret;
    b = bam_init1();
    c = &b->core;
    struct bed_chr *last = NULL;
    int i_bed = -1; // bed iterate
    kstring_t str = {0,0,0};// name buffer
    bam_hdr_t *hdr = args.hdr;
    struct bedaux *bed = args.bed;
    
    while ((ret = sam_read1(args.fp, hdr, b)) >= 0) {
        char *name = hdr->target_name[c->tid];
        if (last == NULL || strcmp(bed->names[last->id], name) != 0) {
            int id = bed_select_chrom(bed, name);
            if (id == -1) {
                last = NULL;
                if (sam_write1(args.out, hdr, b) == -1) error("Failed to write SAM.");
                continue;
            }
            last = &bed->c[id];
            i_bed = 0;
        }
        if (i_bed == -2) { // out range of bed
            if (sam_write1(args.out, hdr, b) == -1) error("Failed to write SAM.");
            continue;
        }
        for (;;) {
            if (i_bed == last->n) break;
            if (last->b[i_bed].end < c->pos+1) i_bed++; // iter bed
            else break;
        }

        if (i_bed == last->n) {
            i_bed = -2;
            if (sam_write1(args.out, hdr, b)== -1) error("Failed to write SAM.");
            continue;
        }
        int end = c->pos + c->l_qseq;
        if (end < last->b[i_bed].start) { // read align before region
            if (sam_write1(args.out, hdr, b) == -1) error("Failed to write SAM.");
            continue;       
        }
        
        uint8_t *tag = bam_aux_get(b, args.tag);
        if (tag) {
            warnings("%s already present at line %s:%d, skip", args.tag, hdr->target_name[c->tid], c->pos+1);
            continue;
        }
        str.l = 0;
        if (last->b[i_bed].name == NULL) {
            ksprintf(&str, "%s:%d-%d", bed->names[last->id], last->b[i_bed].start, last->b[i_bed].end);
        }
        else kputs(last->b[i_bed].name, &str);
        bam_aux_append(b, args.tag, 'Z', str.l+1, (uint8_t*)str.s);
        if (sam_write1(args.out, hdr, b) == -1) error("Failed to write SAM.");
    }
    bam_destroy1(b);
    if (ret != -1) warnings("Truncated file?");
}
/*
void bam_anno_attr_ganno_core(struct ganno *g, bam_hdr_t *h, bam1_t *b)
{
    bam1_core_t *c = &b->core;
    
}
void bam_anno_attr_ganno()
{
    bam1_t *b;
    b = bam_init1();
    int ret;
    bam_hdr_t *hdr = args.hdr;
    
    while ((ret = sam_read1(args.fp, hdr, b)) >= 0) {
        bam_anno_attr_ganno_core(args.ganno, hdr, b);
        if (sam_write1(args.out, hdr, b) == -1) error("Failed to write SAM.");
    }

    bam_destroy1(b);
    if (ret != -1) warnings("Truncated file?");
}
*/
void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.fp);
    sam_close(args.out);
}
int bam_anno_attr(int argc, char *argv[])
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();

    bam_anno_attr_bed();
    //if (args.bed) bam_anno_attr_bed();
    //else bam_anno_attr_ganno();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;    
}


        
