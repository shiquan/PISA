#include "utils.h"
#include "bed.h"
#include "gtf.h"
#include "read_anno.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *gtf_fname;

    int upstream;
    int downstream;

    int stranded;

    struct gtf_spec *G;
    struct bed_spec *B;
    
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .gtf_fname    = NULL,
    .upstream     = 0,
    .downstream   = 0,
    .stranded     = 1,
    .G            = NULL,
    .B            = NULL
};

static int bedanno_usage()
{
    fprintf(stderr, "# Annotate gene element for bed regions.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m annobed -gtf genes.gtf -o anno.bed in.bed\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-gtf  [GTF]     GTF database.\n");
    fprintf(stderr, "-o    [FILE]    Output bed file.\n");
    fprintf(stderr, "-s              Ignore strand.\n");
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

    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-s") == 0) {
            args.stranded = 0;
            continue;
        }

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
    args.G = gtf_read(args.gtf_fname, 2);
    
    return 0;
}

static void summary_report()
{
    
}
static void memory_release()
{
    bed_spec_ext_destroy(args.B);
    bed_spec_destroy(args.B);
    gtf_destroy(args.G);
}

int check_nearest_gene()
{
    return 0;
}

struct anno0 {
    struct gtf *g;
    int type;
    //int dist;
};
// sort by type and distance
static int cmpfunc(const void *_a, const void *_b)
{
    const struct anno0 *a = (const struct anno0*) _a;
    const struct anno0 *b = (const struct anno0*) _b;
    if (a->type != b->type) return (a->type > b->type) - (a->type < b->type);
    return 0;
    // return (a->dist < b->dist) - (a->dist > b->dist);
}

// exon
static int query_exon(int start, int end, struct gtf const *G, struct anno0 *a)
{
    struct anno0 *a0 = malloc(sizeof(struct anno0)* G->n_gtf);
    
    int utr = 0;
    if (G->coding) utr = 1; // if CDS record exists
    
    int i;
    int j = 0;
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf *g0 = G->gtf[i];        
        if (g0->type != feature_exon && g0->type != feature_CDS) continue;
        // debug_print("exon: %d\t%d", g0->start, g0->end);
        if (start > g0->end) continue; // skip this exon
        
        if (start >= g0->start && end <= g0->end) {
            // CDS record is only used to distiguish UTR and EXON
            if (g0->type == feature_CDS) {
                utr = 0;
                continue;
            }
            a0[j].type = BAT_EXON;
            a0[j].g = g0;
            j++;
            continue;
        }
        
        if (end <= g0->start) {
            if (i == 0) break; // first exon

            if (g0->type == feature_CDS) { // before CDS, still treat as UTR..
                continue;
            }
            
            a0[j].type = BAT_INTRON;
            a0[j].g = g0;
            j++;
            break; // no need to check more exons
        }

        // check overlapped records
        
        if (g0->type == feature_CDS) {
            utr = 0; // overlapped with a CDS region; need check type == EXON/EXONINTRON??
            continue;
        }
        
        a0[j].type = BAT_EXONINTRON;
        a0[j].g = g0;
        j++;
    }

    // debug_print("j: %d", j);
    if (j == 0) { // no hit
        a->type = BAT_INTERGENIC;
        a->g = NULL;
        free(a0);
        return 1;
    }

    if (j == 1) {
        a->type = a0[0].type;
        a->g = a0[0].g;
        free(a0);
        return 0;
    }
    
    qsort(a0, j, sizeof(struct anno0), cmpfunc);

    for (i = 0; i < j; ) {
        if (a0[i].type == BAT_INTRON && i > 0) break;
        ++i;
    }

    // debug_print("i: %d", i);
    a->type = a0[0].type;
    a->g =a0[0].g;
    // debug_print("type: %d", a->type);
    // check utr, multiexons
    if (utr == 1) {
        if (G->strand == 0) a->type = BAT_UTR3;
        else a->type = BAT_UTR5;
    } else {
        /* int k; */
        /* for (k = 0; k < i; ++k) { */
        /*     debug_print("%d\t%d\t%d\t%d\t%s", start, end, k, i, bed_typename(a0[k].type)); */
        /* } */
            
        if (i > 1) {
            a->type = BAT_MULTIEXONS;
        }
    }

    free(a0);
    return 0;
}

static int query_trans(int start, int end, struct gtf const *G, struct anno0 *a)
{
    // debug_print("gene: %s, %d\t%d", GTF_genename(args.G,G->gene_name), G->start, G->end);
    struct anno0 *a0 = malloc(sizeof(struct anno0)* G->n_gtf);
    int i;
    int j = 0;
    for (i = 0; i < G->n_gtf; ++i) { // 
        struct gtf *g0 = G->gtf[i];
        if (g0->type != feature_transcript) continue;
        // debug_print("trans: %s, %d\t%d", GTF_transid(args.G,G->gtf[i]->transcript_id), g0->start, g0->end);
        if (start >= g0->end) continue;
        if (end <= g0->start) continue;
        int ret = query_exon(start, end, g0, &a0[j]);
        if (ret == 0) j++;
    }
    
    if (j ==0) { free(a0); return 1; }
    // debug_print("j : %d",j);
    
    if (j > 1) qsort(a0, j, sizeof(struct anno0), cmpfunc);

    a->type = a0[0].type;
    a->g = a0[0].g;
    free(a0);
    // debug_print("type: %d", a->type);
    return 0;
}

int annobed_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bedanno_usage();

    int i;
    for (i = 0; i < args.B->n; ++i) {
        struct bed *b = &args.B->bed[i];
        b->data = NULL;
        
        char *name = dict_name(args.B->seqname, b->seqname);
        struct region_itr *itr = gtf_query(args.G, name, b->start, b->end);
        if (itr == NULL) {
            int id = dict_query(args.G->name, name);
            if (id == -1) error("Chromosome %s not found in GTF, use wrong database?", name);
            
            check_nearest_gene();
            continue;
        }

        // debug_print("itr->n: %d", itr->n);
        // annotate all possibility
        struct anno0 *a = malloc(sizeof(struct anno0)*itr->n);
        // int n = 0;
        int j;
        int k = 0;
        for (j = 0; j < itr->n; ++j) {
            struct gtf const *g0 = (struct gtf*)itr->rets[j];
            if (b->start >= g0->end) continue;
            if (b->end <= g0->start) continue;
            int ret = query_trans(b->start, b->end, g0, &a[k]);
            // debug_print("ret: %d, %d\t%d, bed: %d\t%d", ret, g0->start, g0->end, b->start, b->end);
            if (ret == 0) {
                if (b->strand != -1 && g0->strand != b->strand) {
                    if (a[k].type == BAT_MULTIEXONS) a[k].type = BAT_ANTISENSECOMPLEX;
                    else if (a[k].type == BAT_WHOLEGENE) a[k].type = BAT_ANTISENSECOMPLEX;
                    else if (a[k].type == BAT_UTR3) a[k].type = BAT_ANTISENSEUTR3;
                    else if (a[k].type == BAT_UTR5) a[k].type = BAT_ANTISENSEUTR5;
                    else if (a[k].type == BAT_EXON) a[k].type = BAT_ANTISENSEEXON;
                    else if (a[k].type == BAT_INTRON) a[k].type = BAT_ANTISENSEINTRON;
                    else if (a[k].type == BAT_EXONINTRON) a[k].type = BAT_ANTISENSECOMPLEX;
                }
                k++;
            }
        }

        // debug_print("k, %d",k);
        region_itr_destroy(itr);

        if (k > 1) {
            qsort(a, k, sizeof(struct anno0), cmpfunc);
            int type = a[0].type;
            int j;
            for (j = 0; j < k; ++j) {
                if (a[j].type > type) {
                    k = j;
                    break;
                }
            }
            
        }
        if (k == 0) check_nearest_gene();
        else if (k == 1) {
            // debug_print("%d",a[0].type);                        
            // debug_print("%d\t%d, bed: %d\t%d",a[0].g->start, a[0].g->end, b->start, b->end);

            // debug_print("%s", bed_typename(a[0].type));

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
        } else {
                        
            // debug_print("%d\t%d, bed: %d\t%d",a[0].g->start, a[0].g->end, b->start, b->end);
            // debug_print("%d",a[0].type);
            // debug_print("%s", bed_typename(a[0].type));

            struct bed_ext *e = bed_ext_init();
            e->genes = malloc(k*sizeof(char**));
            e->type = a[0].type;
            int j;
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
            b->data = e;
        }

        free(a);
    }

    bed_spec_write(args.B, args.output_fname, 1);
    
    summary_report();
    memory_release();
    
    return 0;
}
