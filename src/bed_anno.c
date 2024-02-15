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

    int promoter;
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
    .promoter     = 0,
    .upstream     = 2000,
    .downstream   = 100,
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
    fprintf(stderr, "-promoter         Enable promoter regions annotation.\n");
    fprintf(stderr, "-up     [%d]      Define up of TSS as promoter.\n", args.upstream);
    fprintf(stderr, "-down   [%d]      Define downstream of TSS as promoter.\n", args.downstream);
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
        else if (strcmp(a, "-promoter") == 0) {
            args.promoter = 1;
            continue;
        }
        else if (strcmp(a, "-up") == 0) {
            up = a;
            continue;
        }
        else if (strcmp(a, "-down") == 0) {
            down = a;
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

    if (args.promoter) {
        if (up) args.upstream = str2int(up);
        if (down) args.downstream = str2int(down);
    }
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
    int utr = 0;
    if (coding) utr = 1; // if CDS record exists, turn utr to 0 if region overlapped with CDS region
    int pass_cds = 0;

    int n = 0;
    struct gtf **gtf_pool = malloc(G->n_gtf *sizeof(struct gtf*));
    int i;
    for (i = 0; i < G->n_gtf; ++i) { // exon level
        struct gtf *g0 = G->gtf[i];
        // filter type
        if (g0->type != feature_CDS && g0->type != feature_exon) continue;
        // check overlapped
        
        // non-overlap
        if (end <= g0->start) {
            if (n == 0) { // put this record, as intron
                gtf_pool[n++] = G->gtf[i];
            }
            break;   
        }

        if (g0->type == feature_CDS) {
            if (start > g0->end) pass_cds = 1;
        }
            
        if (start > g0->end) continue; // check next
        
        // CDS record is only used to distiguish UTR and EXON
        if (g0->type == feature_CDS) {
            utr = 0; // it's a coding region
            continue;
        }
        
        // push to pool
        gtf_pool[n++] = g0;

        // out of range, no need to check next one
        if (end <= g0->end) break;
    }
    
    if (n == 0) {
        /* if (i > 0 && i != G->n_gtf) { // checked, but no overalpped exons */
        /*     a->type = BAT_INTRON; */
        /* } else { // checked, but out range of transcript */
        a->type = BAT_INTERGENIC;            
        /* } */

        a->g = NULL;
    }
    else if (n == 1) {
        struct gtf *g0 = gtf_pool[0];
        if (g0->start <= start && g0->end >= end) {
            a->type = BAT_EXON;
        }
        else if (end <= g0->start) {
            a->type = BAT_INTRON;
        }
        else {
            a->type = BAT_EXONINTRON;
        }
        a->g = g0;
    }
    else {
        struct gtf *g0 = gtf_pool[0];
        a->type = BAT_MULTIEXONS;
        a->g = g0;
        // struct gtf *g1 = gtf_pool[1];
        // debug_print("%d\t%d\t%d\t%d", g0->start, g0->end, g1->start, g1->end);
    }

    //if ((a->type == BAT_EXON || a->type == BAT_EXONINTRON) && utr == 1) {
    if (a->type == BAT_EXON && utr == 1) {
        // forward
        if (G->strand == 0) a->type = pass_cds ? BAT_UTR3 : BAT_UTR5;
        // backward
        else a->type = pass_cds ? BAT_UTR5 : BAT_UTR3;
    }

    free(gtf_pool);
    // debug_print("%s", bed_typename(a->type));
    
    return a->type == BAT_INTERGENIC;
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
    //debug_print("type: %d", a->type);
    return 0;
}

static int query_promoter(int start, int end, struct gtf *G, struct anno0 *a)
{
    // nonoverlap with this promoter

    int start0 = G->start;
    int end0 = G->start;
    if (G->strand == 1) {
        start0 = G->end;
        end0 = G->end;
        start0 = start0 - args.downstream;
        end0 = start0 + args.upstream;
    } else {
        start0 = start0 - args.upstream;
        end0 = start0 + args.downstream;
    }

    if (start0 < 0) start0 = 1;
    if (end0 < 0) end0 = start0;
    
    if (start >= end0) return 1;
    if (end <= start0) return 1;

    a->type = BAT_PROMOTER;
    a->g = G;
    return 0;
}

int annobed_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bedanno_usage();

    for (int i = 0; i < args.B->n; ++i) {
        struct bed *b = &args.B->bed[i];
        b->data = NULL;
        
        char *name = dict_name(args.B->seqname, b->seqname);
        struct region_itr *itr = gtf_query(args.G, name, b->start, b->end);
        if (itr == NULL) {
            int id = dict_query(args.G->name, name);
            if (id == -1) {
                if (args.skip_chrs == 0) error("Chromosome %s not found in GTF, use wrong database? If not, rerun with `-skip-chrs`.", name);
                b->seqname = -1;
                continue;
            } 
            // check_nearest_gene();
            continue;
        }

        // debug_print("itr->n: %d", itr->n);
        // annotate all possibility
        struct anno0 *a = malloc(sizeof(struct anno0)*itr->n);
        // int n = 0;
        int k = 0;
        for (int j = 0; j < itr->n; ++j) {
            struct gtf *g0 = (struct gtf*)itr->rets[j];
            if (b->start <= g0->start && b->end >= g0->end) {
                a[k].type = BAT_WHOLEGENE;
                a[k].g = g0;
                // k++;
                // continue;
            } else {
                int ret;

                if (args.promoter) {
                    ret = query_promoter(b->start, b->end, g0, &a[k]);
                    if (ret) {
                        ret = query_trans(b->start, b->end, g0, &a[k]);
                    }
                } else {
                    ret = query_trans(b->start, b->end, g0, &a[k]);
                }
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
            //int type = a[0].type;
            for (int j = 0; j < k; ++j) {
                //if (a[j].type > type) {
                if (a[j].type > BAT_WHOLEGENE) {
                    k = j;
                    break;
                }
            }            
        }
        
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
