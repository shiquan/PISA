// Annotate bed files with GTF annotation.
// Required a 6-column bed file, with format as "chr, start, end, name, score, strand".
// Annotated tags such as genes and functional region will be put into extra columns,
// therefore the output file is not exactly a formal bed file; the format of output file is
// "chr, start, end, name, score, strand, n_gene, gene_name, functional region".
// n_gene: the number of gene overlapped
// gene_name: overlapped gene(s), multiple genes seperated by ','
// functional region: unknown, exon, intron, multiexons, exonintron, whole_gene, utr(35),
//                    antisense_utr(35), antisense_intron, antisense_exon, antisense_complex,
//                    intergenic
#include "utils.h"
#include "bed.h"
#include "gtf.h"
#include "read_anno.h"
#include "number.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *gtf_fname;
    const char *report_fname;
    
    int stranded;
    int gene_as_name;
    int skip_chrs;

    /* int promoter; */
    /* int promoter_upstream; */
    /* int promoter_downstream; */

    int upstream;
    int downstream;

    struct gtf_spec *G;
    struct bed_spec *B;
    struct bed_anno_sum *summary;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .gtf_fname    = NULL,
    .report_fname = NULL,
    .stranded     = 1,
    .gene_as_name = 0,
    .skip_chrs    = 0,
    /* .promoter     = 0, */
    /* .promter_upstream     = 2000, */
    /* .promter_downstream   = 100, */
    .upstream  = 1000,
    .downstream= 1000,
    .G            = NULL,
    .B            = NULL,
    .summary      = NULL
};

struct bed_anno_sum {
    uint32_t count;
    uint32_t cov;
};

static int bedanno_usage()
{
    fprintf(stderr, "# Annotate gene element for bed regions.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m annobed -gtf genes.gtf -o anno.bed in.bed\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-gtf    [GTF]     GTF database.\n");
    fprintf(stderr, "-o      [FILE]    Output bed file.\n");
    fprintf(stderr, "-report [FILE]    Summary report. Export in STDERR by default.\n");
    fprintf(stderr, "-is               Ignore strand.\n");
    fprintf(stderr, "-gene-name        Set annatated gene as bed name (column 4).\n");
    fprintf(stderr, "-skip-chrs        Skip chromosomes if not exist in GTF. Defa\n");
    //fprintf(stderr, "-promoter         Enable promoter regions annotation.\n");
    //fprintf(stderr, "-promter-up     [%d]      Define upstream of TSS as promoter.\n", args.promter_upstream);
    //fprintf(stderr, "-promter-down   [%d]      Define downstream of TSS as promoter.\n", args.promter_downstream);
    fprintf(stderr, "-up  [%d]         Annotate intergenic regions at upstream of gene.\n", args.upstream);
    fprintf(stderr, "-down  [%d]       Annotate intergenic regions at downstream of gene.\n", args.downstream);
    fprintf(stderr, "\n\x1b[31m\x1b[1mOutput format\x1b[0m :\n");
    fprintf(stderr, "chromosome,start(0based),end(1based),name,score,strand,number of covered genes, cover gene name(s),type,nearest gene name,distance to nearby gene\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * This tool accepts 3 columns or 6 columns bed file(s), strand (+/-) is set in column 6.\n");
    fprintf(stderr, " * By default, annotation is done with respect to strandness unless -s is set.\n");
    fprintf(stderr, "\n");
    return 1;

}

static int parse_args(int argc, char **argv)
{
    int i;
    const char *up = NULL;
    const char *down = NULL;
    //const char *at_up = NULL;
    //const char *at_down = NULL;
    
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-is") == 0) {
            args.stranded = 0;
            continue;
        }
        else if (strcmp(a, "-gene-name") == 0) {
            args.gene_as_name = 1;
            continue;
        }
        else if (strcmp(a, "-skip-chrs") == 0) {
            args.skip_chrs = 1; 
            continue;
        }
        /* else if (strcmp(a, "-promoter") == 0) { */
        /*     args.promoter = 1; */
        /*     continue; */
        /* } */
        else if (strcmp(a, "-up") == 0) var = &up;
        else if (strcmp(a, "-down") == 0) var = &down;
        /* else if (strcmp(a, "-at-up") == 0) var = &at_up; */
        /* else if (strcmp(a, "-at-down") == 0) var = &at_down; */
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1] != '\0') error("Unknown argument, %s", a);

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument: %s", a);
    }

    if (args.input_fname == NULL) error("No input bed file.");
    if (args.gtf_fname == NULL) error("No gtf file.");

    args.B = bed_read(args.input_fname);
    args.G = gtf_read(args.gtf_fname, 0);

    // if (args.promoter) {
    if (up) args.upstream = str2int(up);
    if (down) args.downstream = str2int(down);
    //}

    /* if (at_up) args.at_upstream = str2int(at_up); */
    /* if (at_down) args.at_downstream = str2int(at_down); */
    
    args.summary = malloc(sizeof(struct bed_anno_sum)*BAT_COUNT);
    memset(args.summary, 0, sizeof(struct bed_anno_sum)*BAT_COUNT);
    
    return 0;
}

static void summary_report()
{
    FILE *out = stderr;
    if (args.report_fname) {
        out = fopen(args.report_fname, "w");
        if (!out) error("%s : %s.", args.report_fname, strerror(errno));
    }
    fputs("Type,Counts,Span(bp)\n", out);
    int i;
    for (i = 1; i < BAT_COUNT; ++i) { // skip unknown
        if (args.summary[i].count > 0)
            fprintf(out, "%s,%d,%d\n", bed_typename(i), args.summary[i].count, args.summary[i].cov);
    }
    fclose(out);
}
static void memory_release()
{
    anno_bed_cleanup();
    free(args.summary);
    bed_spec_ext_destroy(args.B);
    bed_spec_destroy(args.B);
    gtf_destroy(args.G);
}
int annobed_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bedanno_usage();

    for (int i = 0; i < args.B->n; ++i) {
        struct bed *b = &args.B->bed[i];
        b->data = NULL;
        
        char *name = dict_name(args.B->seqname, b->seqname);
        int k;
        struct anno0 *a = anno_bed_core(name, b->start, b->end, b->strand, args.G, &k, args.downstream, args.upstream);
        if (k == 1) {
            struct bed_ext *e = bed_ext_init();
            e->type = a[0].type;
            e->genes = NULL;
            if (a[0].g) {
                e->genes = malloc(sizeof(char**));
                e->genes[0] =  strdup(GTF_genename(args.G, a[0].g->gene_name));
                e->n = 1;
                // debug_print("%s", GTF_genename(args.G, a[0].g->gene_name));
            }
            b->data = e;
        } else if (k > 1) {
            struct bed_ext *e = bed_ext_init();
            struct gtf *g0 = a[0].g;
            int j;
            for (j = 1; j < k; ++j) {
                struct gtf *g1 = a[j].g;
                if (g0->gene_name == g1->gene_name) continue;
                break;
            }

            if (j == k) {
                e->type = a[0].type;
                e->genes = NULL;
                if (a[0].g) {
                    e->genes = malloc(sizeof(char**));
                    e->genes[0] =  strdup(GTF_genename(args.G, a[0].g->gene_name));
                    e->n = 1;
                }
            } else {
                e->genes = malloc(k*sizeof(char**));
                if (a[0].type > 9) {
                    e->type = a[0].type;
                } else {
                    e->type = BAT_MULTIGENES;
                }
                for (j = 0; j < k; ++j) {
                    struct gtf *g = a[j].g;
                    if (g) {
                        e->genes[j] = strdup(GTF_genename(args.G, g->gene_name));
                        // debug_print("%s", GTF_genename(args.G, g->gene_name));
                    } else {
                        e->genes[j]= NULL;
                        error("Should not come here.");
                    }
                }
                e->n = j;
            }
            b->data = e;
        }
        free(a);
    }

    for (int i = 0; i < args.B->n; ++i) {
        struct bed *b = &args.B->bed[i];
        struct bed_ext *e = (struct bed_ext*)b->data;

        if (b->seqname == -1) {
            args.summary[BAT_UNKNOWNCHRS].count++;
            args.summary[BAT_UNKNOWNCHRS].cov+= b->end - b->start;
        } else {
            if (e) {
                args.summary[e->type].count++;
                args.summary[e->type].cov += b->end - b->start;
            } else {
                args.summary[BAT_INTERGENIC].count++;
                args.summary[BAT_INTERGENIC].cov += b->end - b->start;            
            }
        }
    }
    
    bed_spec_write(args.B, args.output_fname, 1, args.gene_as_name);
    
    summary_report();
    memory_release();
    
    return 0;
}
