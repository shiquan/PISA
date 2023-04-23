#include "utils.h"
#include "bed.h"
#include "gtf.h"
#include "read_anno.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *gtf_fname;
    const char *report_fname;
    
    int upstream;
    int downstream;
    
    int stranded;

    int gene_as_name;

    int skip_chrs;
    
    struct gtf_spec *G;
    struct bed_spec *B;
    struct bed_anno_sum *summary;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .gtf_fname    = NULL,
    .report_fname = NULL,
    .upstream     = 0,
    .downstream   = 0,
    .stranded     = 1,
    .gene_as_name = 0,
    .skip_chrs    = 0,
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
    fprintf(stderr, "-s                Ignore strand.\n");
    fprintf(stderr, "-gene-name        Set annatated gene as bed name (column 4).\n");
    fprintf(stderr, "-skip-chrs        Skip chromosomes if not exist in GTF.\n");
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
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-s") == 0) {
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
        fprintf(out, "%s,%d,%d\n", bed_typename(i), args.summary[i].count, args.summary[i].cov);
    }
    fclose(out);
}
static void memory_release()
{
    free(args.summary);
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
static int query_exon(int start, int end, struct gtf const *G, struct anno0 *a, int coding)
{
    struct anno0 *a0 = malloc(sizeof(struct anno0)* G->n_gtf);
    
    int utr = 0;
    if (coding) utr = 1; // if CDS record exists, turn utr==0 if region overlapped with CDS region
    int pass_cds = 0;
    
    int i;
    int j = 0; // iterate for overlapped exons
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf *g0 = G->gtf[i];
        if (g0->type != feature_exon && g0->type != feature_CDS) continue;
        // LOG_print("%s\t%d\t%d\tutr: %d", get_feature_name(g0->type), g0->start, g0->end, utr);
        if (start >= g0->end) {
            // skip CDS
            if (g0->type == feature_CDS) pass_cds = 1; 
            
            continue; // skip this exon
        } 
        
        if (start >= g0->start && end <= g0->end) {
            // CDS record is only used to distiguish UTR and EXON
            if (g0->type == feature_CDS) {
                utr = 0; // it's a coding region
                continue;
            }
            if (j == 0) {
                a0[j].type = BAT_EXON;
                a0[j].g = g0;
            } else if (a0[j-1].type == BAT_EXON) { // not possible ..

                error("Overlapped exon records in one transcript ? %s",  GTF_transid(args.G,G->gtf[i]->transcript_id));
                // a0[j].type = BAT_MULTIEXONS;
                // a0[j].g = g0;
            }
            j++;
            continue; // check if next record is CDS..
        }

        if (j > 0 && a0[j-1].type == BAT_EXON) break; // if already annotated as exon, no need to check downstream regions
        
        if (end <= g0->start) { // overlapped with intron
            if (i == 0) break; // first exon

            if (g0->type == feature_CDS) { // before CDS, still treat as UTR..
                continue; // continue test next record, in case exon record come after CDS
            }
            if (j == 0) { // no exon hit
                a0[j].type = BAT_INTRON;
            }

            // already hit a exon
            else if (a0[j-1].type == BAT_EXON) {
                // should not hit here, if overlapped with more than one exon, last annotation should be exonintron
                error("Should not come here. Type: %s", bed_typename(a0[j-1].type));
                a0[j].type = BAT_MULTIEXONS;
            } else if (a0[j-1].type == BAT_EXONINTRON) {
                a0[j].type = BAT_MULTIEXONS;
            } else if (a0[j-1].type == BAT_MULTIEXONS) {
                
            }
            
            else {
                error("Type: %s", bed_typename(a0[j-1].type));
            }
                
            a0[j].g = g0;
            j++;
            break; // no need to check more exons
        }

        // check overlapped records
        if (g0->type == feature_CDS) {
            utr = 0; // overlapped with a CDS region; need check type == EXON/EXONINTRON??            
            continue; // skip CDS
        }
        
        // region exceed this exon boundary
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

    if (j > 1) qsort(a0, j, sizeof(struct anno0), cmpfunc);

    for (i = 0; i < j; ) {
        if (i > 0 && a0[i].type > BAT_EXONINTRON) break;
        ++i;
    }

    // debug_print("i: %d", i);
    a->type = a0[0].type;
    a->g =a0[0].g;
    // debug_print("type: %d", a->type);
    // check utr, multiexons

    if (a->type == BAT_EXON && utr == 1) {
        // forward
        if (G->strand == 0) a->type = pass_cds ? BAT_UTR3 : BAT_UTR5;
        // backward
        else a->type = pass_cds ? BAT_UTR5 : BAT_UTR3;
    }
    else {
        /* int k; */
        /* for (k = 0; k < i; ++k) { */
        /*     // debug_print("%d\t%d\t%d\t%d\t%s", start, end, k, i, bed_typename(a0[k].type)); */
        /* } */
            
        if (i > 1) { // reset some multi-exonintrons to multiexons
            // LOG_print("Type: %s", bed_typename(a0[0].type));
            a->type = BAT_MULTIEXONS;
        }
    }

    free(a0);
    return 0;
}

static int query_trans(int start, int end, struct gtf const *G, struct anno0 *a)
{
    // nonoverlap with this gene
    if (start >= G->end) return 1;
    if (end <= G->start) return 1;
    
    // debug_print("gene: %s, %d\t%d", GTF_genename(args.G,G->gene_name), G->start, G->end);
    struct anno0 *a0 = malloc(sizeof(struct anno0)* G->n_gtf);
    int i;
    int j = 0;
    for (i = 0; i < G->n_gtf; ++i) { // 
        struct gtf *g0 = G->gtf[i];
        if (g0->type != feature_transcript) continue;
        // debug_print("trans: %s, %d\t%d", GTF_transid(args.G,G->gtf[i]->transcript_id), g0->start, g0->end);

        //  nonoverlap with this trans
        if (start >= g0->end) continue;
        if (end <= g0->start) continue;
        
        int ret = query_exon(start, end, g0, &a0[j], g0->coding);
        if (ret == 0) j++;
        // debug_print("%s", bed_typename(a0[j-1].type));
    }
    
    if (j ==0) { free(a0); return 1; }
    // debug_print("j : %d",j);
    
    if (j > 1) qsort(a0, j, sizeof(struct anno0), cmpfunc);

    // if hit more than one exon/cds record, assign by the first record
    a->type = a0[0].type;
    a->g = a0[0].g;
    // debug_print("%s\t%s", bed_typename(a->type), GTF_transid(args.G, a->g->transcript_id));
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
            if (id == -1) {
                if (args.skip_chrs == 0) error("Chromosome %s not found in GTF, use wrong database?", name);
                b->seqname = -1;
                continue;
            } 
            
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
            struct gtf *g0 = (struct gtf*)itr->rets[j];
            if (b->start <= g0->start && b->end >= g0->end) {
                a[k].type = BAT_WHOLEGENE;
                a[k].g = g0;
                // k++;
                // continue;
            } else {            
                int ret = query_trans(b->start, b->end, g0, &a[k]);
                // debug_print("ret: %d, %d\t%d, bed: %d\t%d", ret, g0->start, g0->end, b->start, b->end);
                if (ret != 0) continue; // nonoverlapped
            }
            
            // for stranded, reannotate antisense 
            if (b->strand != -1 && g0->strand != b->strand) {
                if (a[k].type == BAT_MULTIEXONS) a[k].type = BAT_ANTISENSECOMPLEX;
                else if (a[k].type == BAT_WHOLEGENE) a[k].type = BAT_ANTISENSECOMPLEX;
                else if (a[k].type == BAT_EXONINTRON) a[k].type = BAT_ANTISENSECOMPLEX;
                else if (a[k].type == BAT_UTR3) a[k].type = BAT_ANTISENSEUTR3;
                else if (a[k].type == BAT_UTR5) a[k].type = BAT_ANTISENSEUTR5;
                else if (a[k].type == BAT_EXON) a[k].type = BAT_ANTISENSEEXON;
                else if (a[k].type == BAT_INTRON) a[k].type = BAT_ANTISENSEINTRON;
            }
            k++;
            
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
            e->type = BAT_MULTIGENES;
            
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

        struct bed_ext *e = (struct bed_ext*)b->data;
        
        if (e) {
            args.summary[e->type].count++;
            args.summary[e->type].cov += b->end - b->start;
        } else {
            args.summary[BAT_INTERGENIC].count++;
            args.summary[BAT_INTERGENIC].cov += b->end - b->start;            
        }
        free(a);
    }

    bed_spec_write(args.B, args.output_fname, 1, args.gene_as_name);
    
    summary_report();
    memory_release();
    
    return 0;
}
