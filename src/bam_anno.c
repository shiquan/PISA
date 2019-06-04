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
#include "gtf.h"
#include <zlib.h>

static int usage()
{
    fprintf(stderr, "* Annotate bam records with overlapped function regions. Such as gene, trnascript etc.\n");
    fprintf(stderr, "anno_bam  -bed peak.bed -tag PK -o anno.bam in.bam\n");
    fprintf(stderr, "Options :\n");
    fprintf(stderr, "  -o               Output bam file.\n");
    fprintf(stderr, "Options for BED file :\n");
    fprintf(stderr, "  -bed             Function regions. Three or four columns bed file. Col 4 could be empty or names of this region.\n");
    fprintf(stderr, "  -tag             Attribute tag name. Set with -bed\n");
    fprintf(stderr, "Options for GTF file :\n");
    fprintf(stderr, "  -gtf             GTF annotation file. -gtf is conflict with -bed, if set strand will be consider.\n");
    fprintf(stderr, "  -tags            Attribute names. Default is TX,AN,GN,GX,RE.\n");
    fprintf(stderr, "  -ignore-strand   Ignore strand of transcript in GTF. Reads mapped to antisense transcripts will also be count.\n");
    fprintf(stderr, "Notice : * For GTF mode, this program will set tags in default, you could also reset them by -tags.\n");
    fprintf(stderr, "           TX : Transcript id.\n");
    fprintf(stderr, "           AN : Same with TX but set only if read mapped to antisense strand of transcript.\n");
    fprintf(stderr, "           GN : Gene name.\n");
    fprintf(stderr, "           GX : Gene ID.\n");
    fprintf(stderr, "           RE : Region type, should E(exonic), N(intronic), I(intergenic)\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bed_fname;
    const char *tag; // attribute in BAM
    const char *gtf_fname;
    //    const char *tags;
    // int n_tag;
    char **tags;
    int ignore_strand;

    htsFile *fp;
    htsFile *out;
    bam_hdr_t *hdr;

    struct gtf_spec *G;

    struct bedaux *B;
    struct bed_chr *last;
    int i_bed;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .bed_fname = NULL,
    
    .tag = NULL,
    
    .gtf_fname = NULL,
    .tags =NULL,
    .ignore_strand = 0,

    .fp = NULL,
    .out = NULL,
    .hdr = NULL,

    .G = NULL,
    .B = NULL,
    .last = NULL,
    .i_bed = -1,
};
static int parse_args(int argc, char **argv)
{
    int i;
    const char *tags = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-tags") == 0) var = &tags;
        else if (strcmp(a, "-ignore-strand") == 0) {
            args.ignore_strand = 1;
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
        error("Unknown argument: %s", a);
    }
    // CHECK_EMPTY(args.bed_fname, "-bed must be set.");

    if (args.bed_fname == NULL && args.gtf_fname == NULL) 
        error("-bed or -gtf must be set.");
    if (args.bed_fname && args.gtf_fname)
        error("-bed is conflict with -gtf, you can only choose one mode.");
    if (args.ignore_strand && args.gtf_fname == NULL)
        error("Only set -ignore-strand with -gtf.");
    
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.tag, "-tag must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");

    args.fp  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    args.hdr = sam_hdr_read(args.fp);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    //int n_bed = 0;
    
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    if (args.bed_fname) {
        args.B = bed_read(args.bed_fname);
        if (args.B == 0 || args.B->n == 0) error("Bed is empty.");
    }
    else {
        args.G = gtf_read(args.gtf_fname);
        if (args.G == NULL) error("GTF is empty.");
    }
    
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");

    return 0;
}

int check_is_overlapped_bed(bam_hdr_t *hdr, bam1_t *b, struct bedaux *B)
{
    bam1_core_t *c;
    c = &b->core;
    char *name = hdr->target_name[c->tid];
    if (args.last == NULL || strcmp(B->names[args.last->id], name) != 0) {
        int id = bed_select_chrom(B, name);
        if (id == -1) {
            args.last = NULL;
            return 0;
        }

        args.last = &B->c[id];
        args.i_bed = 0;
    }
    if (args.i_bed == -2) { // out range of bed
        return 0;
    }
    
    for (;;) {
        if (args.i_bed == args.last->n) break;
        if (args.last->b[args.i_bed].end < c->pos+1) args.i_bed++; // iter bed
        else break;
    }

    if (args.i_bed == args.last->n) {
        args.i_bed = -2;
        return 0;
    }
    int end = c->pos + c->l_qseq;
    if (end < args.last->b[args.i_bed].start) { // read align before region
        return 0;
    }
    
    uint8_t *tag = bam_aux_get(b, args.tag);
    if (tag) {
        warnings("%s already present at line %s:%d, skip", args.tag, hdr->target_name[c->tid], c->pos+1);
        return 1;
    }
    kstring_t str = {0,0,0};// name buffer
    if (args.last->b[args.i_bed].name == NULL) {
        ksprintf(&str, "%s:%d-%d", B->names[args.last->id], args.last->b[args.i_bed].start, args.last->b[args.i_bed].end);
    }
    else kputs(args.last->b[args.i_bed].name, &str);
    bam_aux_append(b, args.tag, 'Z', str.l+1, (uint8_t*)str.s);
    free(str.s);
    return 0;
}

int check_is_overlapped_gtf(bam_hdr_t *h, bam1_t *b, struct gtf_spec *G)
{
    bam1_core_t *c;
    c = &b->core;
    char *name = h->target_name[c->tid];
    struct gtf_lite *g;
    int n = 0;
    // todo: improve the algorithm of overlap
    g = gtf_overlap(G, name, c->pos, c->pos+c->l_qseq, &n);
    if (n==0) return 1;
    if (c->flag & BAM_FREVERSE) {
        if (g->strand == 0) goto anitisense;        
    }
    else if (g->strand == 1) goto antisense;

    int l;
    // TX
    char *TX_tag = get_TX_tag();    
    char *trans = gtf_get_transcript(g);
    l = strlen(trans);
    bam_aux_append(b, AN_tag, 'Z', l+1, (uint8_t*)trans);

    // GN
    char *GN_tag = get_GN_tag();
    char *gene = gtf_get_gene(g);
    l = strlen(gene);
    bam_aux_append(b, GN_tag, 'Z', l+1, (uint8_t*)gene);

    // GX
    char *GX_tag = get_GX_tag();
    char *gene_id = gtf_get_geneid(g);
    l = strlen(gene_id);
    bam_aux_append(b, GX_tag, 'Z', l+1, (uint8_t*)gene_id);
    // RE
    char *RE_tag = get_RE_tag();
    int type = get_get_type(g, b);
    bam_aux_append(b, RE_tag, 'A', 1, type==0? 'E' : 'I');
    
    return 0;
    

  antisense:
    // AN
    char *AN_tag = get_AN_tag();
    char *trans = gtf_get_transcript(g);
    bam_aux_append(b, AN_tag, 'Z', str.l+1, (uint8_t*)trans);
    return 0;
}
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
    
    bam1_t *b;
    int ret;
    b = bam_init1();
    
    while ((ret = sam_read1(args.fp, args.hdr, b)) >= 0) {

        if (args.B) 
            check_is_overlapped_bed(args.hdr, b, args.B); 
        else
            check_is_overlapped_gtf(args.hdr, b, args.G);
        
        if (sam_write1(args.out, args.hdr, b) == -1) error("Failed to write SAM.");
    }


    bam_destroy1(b);

    if (ret != -1) warnings("Truncated file?");

    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;    
}


        
