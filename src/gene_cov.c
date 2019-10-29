#include "utils.h"
#include "dict.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "number.h"
#include "gtf.h"
#include "htslib/kstring.h"

static int usage()
{
    fprintf(stderr, "* Stat reads coverage over gene body for each cell.\n");
    fprintf(stderr, "gene_cov [options] in.bam\n");
    fprintf(stderr, "  -gtf         GTF file.\n");
    fprintf(stderr, "  -unit        Stat unit, default is gene. [gene|transcript|exon|cds]\n");
    fprintf(stderr, "  -tag         Cell barcode tag.\n");
    fprintf(stderr, "  -list        Cell barcode white list.\n");
    fprintf(stderr, "  -gene        Gene white list.\n");
    fprintf(stderr, "  -trans       Transcript white list.\n");
    fprintf(stderr, "  -summary     Summary output.\n");
    fprintf(stderr, "  -bulk        Bulk coverage table.\n");
    fprintf(stderr, "  -o           Cell barcode X Gene coverage matrix.\n");
    fprintf(stderr, "  -@           Thread to unpack BAM.\n");
    return 1;
}

void *debug_malloc(int x)
{
    debug_print("Try to allocate %d bits.", x);
    return malloc(x);
}
void *debug_realloc(void *a, int x)
{
    debug_print("Try to re-allocate %d bits.", x);
    return realloc(a, x);
}
// Output
// Cell_barcode Gene Coverage

// Summary output
// Coverge  Gene_counts Cell_barcode

static struct args {
    const char *input_fname;
    const char *gtf_fname;
    const char *cell_list;
    const char *gene_list;
    const char *transcript_list;
    const char *summary_fname;
    const char *output_fname;
    const char *bulk_fname;
    const char *tag;
    int file_th;

    htsFile    *fp;
    hts_idx_t  *idx;
    bam_hdr_t  *hdr;
    FILE       *fp_bulk;
    FILE       *fp_mtx;
    FILE       *fp_summary;
    struct gtf_spec *G;
} args = {
    .input_fname     = NULL,
    .gtf_fname       = NULL,
    .cell_list       = NULL,
    .gene_list       = NULL,
    .transcript_list = NULL,
    .summary_fname   = NULL,
    .output_fname    = NULL,
    .bulk_fname      = NULL,
    .file_th         = 4,
    .tag             = NULL,
    .fp              = NULL,
    .idx             = NULL,
    .hdr             = NULL,
    .fp_bulk         = NULL,
    .fp_mtx          = NULL,
    .fp_summary      = NULL,
    .G               = NULL,
};

enum stat_unit {
    unit_gene,
    unit_transcript,
    unit_exon,
    unit_cds,
};

static enum stat_unit stat_unit = unit_gene;

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *file_th = NULL;
    const char *unit = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-list") == 0) var = &args.cell_list;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-gene") == 0) var = &args.gene_list;
        else if (strcmp(a, "-summary") == 0) var = &args.summary_fname;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-bulk") == 0) var = &args.bulk_fname;
        else if (strcmp(a, "-unit") == 0) var = &unit;
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
    if (args.gtf_fname == NULL) error("-gtf must be set.");
    if (args.cell_list == NULL) error("Cell barcode list should be set with -list.");

    if (file_th) args.file_th = str2int(file_th);

    if (unit) {
        if (strcmp(unit, "gene") == 0) stat_unit = unit_gene;
        else if (strcmp(unit, "transcript") == 0) stat_unit = unit_transcript;
        else if (strcmp(unit, "exon") == 0) stat_unit = unit_exon;
        else if (strcmp(unit, "cds") == 0) stat_unit = unit_cds;
        else error("Unknown stat unit, should be [gene|transcript|exon|cds]. %s", unit);
    }

    
    args.fp  = hts_open(args.input_fname, "r");
    if (args.fp == NULL)
        error("%s : %s.", args.input_fname, strerror(errno));
    
    htsFormat type = *hts_get_format(args.fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    args.idx = sam_index_load(args.fp, args.input_fname);
    if (args.idx == NULL)
        error("BAM file need be indexed.");

    args.hdr = sam_hdr_read(args.fp);
    CHECK_EMPTY(args.hdr, "Failed to open header.");

    hts_set_threads(args.fp, args.file_th);

    args.fp_bulk = args.bulk_fname == NULL ? stdout : fopen(args.bulk_fname, "w");
    if (args.fp_bulk == NULL) error("%s : %s.", args.bulk_fname, strerror(errno));
    
    args.G = gtf_read(args.gtf_fname, 1);
    if (args.G == NULL) error("GTF is empty.");

    if (args.output_fname) {
        args.fp_mtx = fopen(args.output_fname, "w");
        if (args.fp_mtx == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    }

    if (args.summary_fname) {
        args.fp_summary = fopen(args.summary_fname, "w");
        if (args.fp_summary == NULL) error("%s : %s.", args.summary_fname, strerror(errno));
    }
    return 0;
}

static void memory_release()
{
    fclose(args.fp_bulk);
    if (args.fp_mtx) fclose(args.fp_mtx);
    if (args.fp_summary) fclose(args.fp_summary);
    hts_idx_destroy(args.idx);
    bam_hdr_destroy(args.hdr);
    gtf_destory(args.G);
    hts_close(args.fp);
}

struct bed {
    int start;
    int end;
};
struct cov {
    int n, m;
    struct bed *bed;
};

struct gcov {
    struct dict *names;
    struct dict *bcodes; // load from barcode list    
    struct cov *temp_cov; // temp coverage record for last gene
    // uint8_t **cov;
};

static int cmpfunc(const void *_a, const void *_b)
{
    struct bed *a = (struct bed *)_a;
    struct bed *b = (struct bed *)_b;
    
   return a->start - b->start == 0 ? a->end - b->end : a->start - b->start;
}
static void cov_merge1(struct cov *cov)
{
    assert(cov->n > 0);
    qsort(cov->bed, cov->n, sizeof(struct bed), cmpfunc);
    int st = cov->bed[0].start;
    int ed = cov->bed[0].end;
    int i;
    int j = 0;
    for (i = 1; i < cov->n; ++i) {
        struct bed *bed = &cov->bed[i];
        assert(st <= bed->start);
        if (ed >= bed->start) {
            if (ed < bed->end) ed = bed->end;
        }
        else {
            cov->bed[j].start = st;
            cov->bed[j].end = ed;
            j++;
            st = bed->start;
            ed = bed->end;
        }
    }
    cov->bed[j].start = st;
    cov->bed[j].end = ed;
    ++j;
    cov->n = j;
}
int cov_sum(struct cov *cov)
{
    if (cov->n == 0) return 0;

    cov_merge1(cov);
    int sum = 0;
    int i;
    for (i = 0; i < cov->n; ++i) sum += cov->bed[i].end - cov->bed[i].start + 1;
    return sum;
}
int cov_sum2(struct cov *cov1, struct cov *cov2)
{
    struct cov *cov = malloc(sizeof(struct cov));
    cov->m = cov1->n + cov2->n;
    cov->n = 0;
    //kroundup32(cov->m);
    cov->bed = malloc(cov->m *sizeof(struct bed));
    int i;
    for (i = 0; i < cov1->n; ++i) {
        cov->bed[i+cov->n].start = cov1->bed[i].start;
        cov->bed[i+cov->n].end   = cov1->bed[i].end;
    }
    cov->n = cov1->n;
    for (i = 0; i < cov2->n; ++i) {
        cov->bed[i+cov->n].start = cov2->bed[i].start;
        cov->bed[i+cov->n].end   = cov2->bed[i].end;
    }
    cov->n += cov2->n;
    int sum;
    sum = cov_sum(cov);
    free(cov->bed);
    free(cov);
    return sum;
}

struct gcov *gcov_init(const char *fn)
{
    struct gcov *g = malloc(sizeof(*g));
    memset(g, 0, sizeof(*g));
    g->names = dict_init();
    g->bcodes = dict_init();
    if (dict_read(g->bcodes, fn))
        error("Failed to load cell barcode.");

    g->temp_cov = malloc(dict_size(g->bcodes)*sizeof(struct cov));
    memset(g->temp_cov, 0, dict_size(g->bcodes)*sizeof(struct cov));
    /*
    g->cov = malloc(dict_size(g->bcodes)*sizeof(void*));
    int i;
    for (i = 0; i < dict_size(g->bcodes); ++i) {
        g->cov[i] = malloc(101);
        memset(g->cov[i], 0, 101);
    }
    */
    return g;
}
void gcov_destory(struct gcov *g)
{
    int i;
    // for (i = 0; i < dict_size(g->bcodes); ++i) free(g->cov[i]);
    for (i = 0; i < dict_size(g->bcodes); ++i) free(g->temp_cov[i].bed);
    free(g->temp_cov);
    // free(g->cov);
    dict_destroy(g->names);
    dict_destroy(g->bcodes);
    free(g);
}
static void cov_push(struct cov *cov, int start, int end)
{
    if (cov->n == cov->m) {
        if (cov->n > 1000) cov_merge1(cov);
        if (cov->n == cov->m) {
            cov->m = cov->m == 0 ? 100 : cov->m<<1;
            cov->bed = realloc(cov->bed, cov->m*sizeof(struct bed));
        }
    }
    cov->bed[cov->n].start = start;
    cov->bed[cov->n].end = end;
    cov->n++;
}

static void gl2bed(struct cov *cov, struct gtf_lite *gl)
{
    if (gl->n_son > 0) {
        int i;
        for (i = 0; i < gl->n_son; ++i) {
            struct gtf_lite *g1 = &gl->son[i];
            gl2bed(cov,g1);
        }
    }
    else {
        cov_push(cov, gl->start, gl->end);
    }    
}

struct acc_gene_cov {
    int n, m;
    uint8_t *cov;
};
int write_summary(struct gcov *g,struct acc_gene_cov *acc_gene_cov)
{
    if (args.fp_summary == NULL) return 1;
    if (acc_gene_cov->n == 0) return 1;
    
    int i;
    int count[101];
    memset(count, 0, 101*sizeof(int));
    for (i = 0; i < acc_gene_cov->n; ++i) count[acc_gene_cov->cov[i]]++;
    for (i = 0; i < 101; ++i)  
        fprintf(args.fp_summary,"Accumulation\t%d\t%d\n", i, count[i]);
    /*
    for (i = 0; i < dict_size(g->bcodes); ++i) {
        int j;
        for (j = 0; j < 101; ++j)
            fprintf(args.fp_summary,"%s\t%d\t%d\n", dict_name(g->bcodes, i), j, g->cov[i][j]);
    }
    */
    return 0;
}
int gene_cov_core(htsFile *fp, hts_idx_t *idx, char *name, int tid, struct gtf_lite *gl, struct gcov *gcov, struct acc_gene_cov *acc_gene_cov)
{
    bam1_t *b = bam_init1();
    struct cov *gene_bed = malloc(sizeof(struct cov));
    memset(gene_bed, 0, sizeof(struct cov));
    gl2bed(gene_bed, gl);
    int lgen = cov_sum(gene_bed);
    
    struct cov *acc_cov = malloc(sizeof(struct cov));
    memset(acc_cov, 0, sizeof(struct cov));
    // reset buffer
    int k;
    for (k = 0; k < dict_size(gcov->bcodes); ++k) gcov->temp_cov[k].n = 0;
    
    int ex;
    for (ex = 0; ex < gene_bed->n; ++ex) {
        //hts_itr_t *itr = sam_itr_queryi(idx, tid, gl->start, gl->end);
        struct bed *exon = &gene_bed->bed[ex];
        hts_itr_t *itr = sam_itr_queryi(idx, tid, exon->start, exon->end);
        int r;
        
        // each block come from exactly one gene
        while ((r = sam_itr_next(fp, itr, b)) >= 0) {
            uint8_t *tag = bam_aux_get(b, args.tag);
            if (!tag) continue;
            char *cb = (char*)(tag+1);
            int id;
            id = dict_query(gcov->bcodes, cb);
            if (id == -1) continue;

            // for counting accumulation coverage
            int l = 0;
            int start = b->core.pos+1;        
            dict_push(gcov->names, name);
            
            struct cov *cov = &gcov->temp_cov[id];
            int k;
            for (k = 0; k < b->core.n_cigar; ++k) {
                int cig = bam_cigar_op(bam_get_cigar(b)[k]);
                int ncig = bam_cigar_oplen(bam_get_cigar(b)[k]);
                if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
                    l += ncig;
                }
                else if (cig == BAM_CDEL) {
                    l += ncig;
                }
                else if (cig == BAM_CREF_SKIP) {
                    //cov_push(cov, start, start+l-1);
                    //cov_push(acc_cov, start, start+l-1);
                    // reset block
                    start = start + l + ncig;
                    l = 0;
                }
            }
            if (l != 0) {
                cov_push(cov, start, start+l-1);
                cov_push(acc_cov, start, start+l-1);
            }
       
        }
         
        hts_itr_destroy(itr);
    }
    
    int lcov = cov_sum(acc_cov);
    if (lcov == 0) goto not_update_cov;
    
    // count coverage of this block and reset buffer
    int gid = dict_query(gcov->names, name);
    if (gid == -1) goto not_update_cov;
    if (gid != dict_size(gcov->names)-1) goto not_update_cov;
    
    kstring_t str = {0,0,0};
    kputs(name, &str);    
    for (k = 0; k < dict_size(gcov->bcodes); ++k) {
        struct cov *cov = &gcov->temp_cov[k];        
        int lcov = cov_sum(cov);
        uint8_t r = 0;
        if (lcov != 0) {
            int sum = cov_sum2(gene_bed, cov);
            float f = (float)(lgen + lcov - sum)/lgen;
            r = (uint8_t)(f*100);
        }
        // gcov->cov[k][r]++;
        kputc('\t', &str);
        kputw(r, &str);
    }

    kputc('\n', &str);
    if (args.fp_mtx) fputs(str.s, args.fp_mtx);
    free(str.s);
   
    // for (k = 0; k < dict_size(gcov->bcodes); ++k) gcov->temp_cov[k].n = 0;

    /*
    if (acc_gene_cov->n == acc_gene_cov->m) {
        acc_gene_cov->m = acc_gene_cov->m == 0 ? 1024 : acc_gene_cov->m << 1;
        acc_gene_cov->cov = realloc(acc_gene_cov->cov, acc_gene_cov->m);            
    }
    */
    int sum  = cov_sum2(gene_bed, acc_cov);
    assert(lgen+lcov-sum <= lgen);
    // acc_gene_cov->cov[acc_gene_cov->n] = (uint8_t)(((float)(lgen+lcov-sum)/lgen)*100);
    
    //fprintf(args.fp_bulk, "%s\t%d\t%d\n", name, lgen, acc_gene_cov->cov[acc_gene_cov->n]);
    //acc_gene_cov->n++;
    fprintf(args.fp_bulk, "%s\t%d\t%d\n", name, lgen, (uint8_t)(((float)(lgen+lcov-sum)/lgen)*100));
    free(acc_cov->bed);
    free(acc_cov);
    free(gene_bed->bed);
    free(gene_bed);
    bam_destroy1(b);
    return 0;

  not_update_cov:
    bam_destroy1(b);
    free(gene_bed->bed);
    free(gene_bed);
    free(acc_cov->bed);
    free(acc_cov);
    return 1;
}
int gene_cov(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    struct gcov *gcov = gcov_init(args.cell_list);

    struct dict *gene_dict = NULL;
    if (args.gene_list) {
        gene_dict = dict_init();
        if (dict_read(gene_dict, args.gene_list))
            error("Failed to load gene list.");
    }

    struct acc_gene_cov *acc_gene_cov;
    acc_gene_cov = malloc(sizeof(struct acc_gene_cov));
    memset(acc_gene_cov, 0, sizeof(struct acc_gene_cov));

    struct gtf_spec *G = args.G;

    int i;
    kstring_t str = {0,0,0};
    kputs("Name", &str);
    for (i = 0; i < dict_size(gcov->bcodes); ++i) {
        kputc('\t', &str);
        kputs(dict_name(gcov->bcodes, i), &str);
    }

    if (args.fp_mtx) {
        fputs(str.s, args.fp_mtx);
        fputc('\n', args.fp_mtx);
    }

    free(str.s);
    
    // iter from whole gtf
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *gl = &G->gtf[i];
        if (gene_dict) {
            int ret = dict_query(gene_dict, dict_name(G->gene_name, gl->gene_name));
            if (ret == -1) continue;
        }
        
        int tid;       
        tid = bam_name2id(args.hdr,dict_name(G->name, gl->seqname));
        if (tid == -1) continue;

        char *gene = dict_name(G->gene_name, gl->gene_name);        

        if (stat_unit == unit_gene) {
            gene_cov_core(args.fp, args.idx, gene, tid, gl, gcov, acc_gene_cov);
        }
        else if (stat_unit == unit_transcript) {
            int j;
            for (j = 0; j < gl->n_son; ++j) {
                struct gtf_lite *g1 = &gl->son[j];
                if (g1->type != feature_transcript) continue;
                char *name = dict_name(G->transcript_id, g1->transcript_id);
                assert (name != NULL);
                gene_cov_core(args.fp, args.idx, name, tid, g1, gcov, acc_gene_cov);
            }
        }
        else if (stat_unit == unit_exon) {
            int j;
            for (j = 0; j < gl->n_son; ++j) {
                struct gtf_lite *g1 = &gl->son[j];
                if (g1->type != feature_transcript) continue;
                int k;
                for (k = 0; k < g1->n_son; ++k) {
                    struct gtf_lite *g2 = &g1->son[k];
                    if (g2->type != feature_exon) continue;
                    kstring_t str = {0,0,0};
                    kputs(dict_name(G->name,g2->seqname), &str);
                    kputc(':', &str);
                    kputw(g2->start, &str);
                    kputc('-', &str);
                    kputw(g2->end, &str);
                    gene_cov_core(args.fp, args.idx, str.s, tid, g2, gcov, acc_gene_cov);
                    free(str.s);
                }
            }
        }
        else if (stat_unit == unit_cds) {
            int j;
            for (j = 0; j < gl->n_son; ++j) {
                struct gtf_lite *g1 = &gl->son[j];
                if (g1->type != feature_transcript) continue;
                int k;
                for (k = 0; k < g1->n_son; ++k) {
                    struct gtf_lite *g2 = &g1->son[k];
                    if (g2->type != feature_CDS) continue;
                    kstring_t str = {0,0,0};
                    kputs(dict_name(G->name,g2->seqname), &str);
                    kputc(':', &str);
                    kputw(g2->start, &str);
                    kputc('-', &str);
                    kputw(g2->end, &str);
                    gene_cov_core(args.fp, args.idx, str.s, tid, g2, gcov, acc_gene_cov);
                    free(str.s);
                }
            }            
        }
    }

    // write_cov_mtx(args.output_fname, gcov);
    write_summary(gcov, acc_gene_cov);

    if (acc_gene_cov->m) free(acc_gene_cov->cov);
    free(acc_gene_cov);
    gcov_destory(gcov);
    if (gene_dict) dict_destroy(gene_dict);
    memory_release();
    return 0;
}
