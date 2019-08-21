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
    fprintf(stderr, "  -tag         Cell barcode tag.\n");
    fprintf(stderr, "  -list        Cell barcode white list.\n");
    fprintf(stderr, "  -gene        Gene white list.\n");
    fprintf(stderr, "  -trans       Transcript white list.\n");
    fprintf(stderr, "  -summary     Summary output.\n");
    fprintf(stderr, "  -o           Cell barcode X Gene coverage matrix.\n");
    fprintf(stderr, "  -@           Thread to unpack BAM.\n");
    return 1;
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

    const char *tag;
    int file_th;
} args = {
    .input_fname     = NULL,
    .gtf_fname       = NULL,
    .cell_list       = NULL,
    .gene_list       = NULL,
    .transcript_list = NULL,
    .summary_fname   = NULL,
    .output_fname    = NULL,
    .file_th         = 4,
    .tag             = NULL,
};

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *file_th = NULL;

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
    
    return 0;
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
    struct dict *genes;
    struct dict *bcodes; // load from barcode list    
    struct cov *bed; // temp coverage record for last gene
    int n, m;
    uint8_t **covs;
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
    kroundup32(cov->m);
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

    cov_merge1(cov);
    int sum;
    sum = cov_sum(cov);
    free(cov->bed);
    free(cov);
    return sum;
}

struct gcov *gcov_init()
{
    struct gcov *g = malloc(sizeof(*g));
    memset(g, 0, sizeof(*g));
    g->genes = dict_init();
    g->bcodes = dict_init();
    return g;
}
void gcov_destory(struct gcov *g)
{
    int i;
    for (i = 0; i < dict_size(g->genes); ++i) free(g->covs[i]);
    free(g->covs);
    dict_destroy(g->genes);
    dict_destroy(g->bcodes);
    free(g);
}
static void cov_push(struct cov *cov, int start, int end)
{
    if (cov->n == cov->m) {
        cov->m = cov->m == 0 ? 1024 : cov->m <<1;
        kroundup32(cov->m);
        cov->bed = realloc(cov->bed, cov->m*sizeof(struct bed));
    }
    cov->bed[cov->n].start = start;
    cov->bed[cov->n].end = end;
    cov->n++;    
}
static int push_gene_cov(struct gcov *g, char *gene, char *CB, bam1_t *b)
{
    int ret;
    ret = dict_query(g->genes, gene);
    if (ret == -1) {
        ret = dict_push(g->genes, gene);
        g->bed = malloc(dict_size(g->genes)*sizeof(struct cov));
        memset(g->bed, 0, dict_size(g->genes)*sizeof(struct cov));
        if (g->n == g->m) {
            g->m = g->m == 0 ? 1024 : g->m<<1;
            g->covs = realloc(g->covs, g->m*sizeof(void*));
        }

        if (g->n < dict_size(g->genes)) {
            int i;
            for (i = g->n; i < dict_size(g->genes); ++i) {
                g->covs[i] = malloc(dict_size(g->bcodes)*1);
                memset(g->covs[i], 0, dict_size(g->bcodes));
            }
            g->n = dict_size(g->genes);
        }
    }

    int id = dict_query(g->bcodes, CB);
    struct cov *bed = &g->bed[id];
    
    int endpos = bam_endpos(b);
    cov_push(bed, b->core.pos+1, endpos);

    return 0;
}
static struct cov *gl2bed(struct gtf_lite *gl)
{
    struct cov *cov = malloc(sizeof(struct cov));
    memset(cov, 0, sizeof(struct cov));
    int i;
    for (i = 0; i < gl->n_son; ++i) {
        struct gtf_lite *g1 = &gl->son[i];
        assert(g1->type == feature_transcript);
        int j;
        for (j = 0; j < g1->n_son; ++j) {
            struct gtf_lite *g2 = &g1->son[i];
            cov_push(cov, g2->start, g2->end);
        }
    }
    return cov;
}
int calc_buf_gene_cov(struct gcov *g, struct gtf_spec *G, struct gtf_lite *gl)
{
    int gid = dict_query(g->genes, dict_name(G->gene_name,gl->gene_name));
    assert(gid>=0);
    struct cov *gene_bed = gl2bed(gl);    
    int i;
    int lgen = cov_sum(gene_bed);
    // debug_print("lgene : %d", lgen);        
    for (i = 0; i < dict_size(g->bcodes); ++i) {
        g->covs[gid][i] = 0;
        struct cov *cov = &g->bed[gid];
        int lcov = cov_sum(cov);
        // debug_print("%s : %d", dict_name(g->bcodes, i), lcov);        
        if (lcov == 0) continue;
        int sum = cov_sum2(gene_bed, cov);
        float f = (float)(lgen + lcov - sum)/lgen;
        g->covs[gid][i] = (uint8_t)(f*100);        
    }
    return 0;
}
int write_cov_mtx(const char *fname, struct gcov *g)
{
    if (fname == NULL) return 1;
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
        error("%s : %s.", fname, strerror(errno));

    int i;
    fputs("Gene", fp);
    for (i = 0; i < dict_size(g->bcodes); ++i) {
        fputc('\t', fp);
        fputs(dict_name(g->bcodes, i), fp);
    }
    fputc('\n', fp);

    int j;
    for (j = 0; j < dict_size(g->genes); ++j) {
        fputs(dict_name(g->genes, j), fp);
        for (i = 0; i < dict_size(g->bcodes); ++i)
            fprintf(fp,"\t%d", g->covs[j][i]);

        fputc('\n', fp);
    }
    fclose(fp);
    return 0;
}
int write_summary(const char *fname, struct gcov *g)
{
    if (fname == NULL) return 1;
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
        error("%s : %s.", fname, strerror(errno));
    
    int i;
    for (i = 0; i < dict_size(g->bcodes); ++i) {
        int count[101];
        memset(count, 0, 101*sizeof(int));
        int j;
        for (j = 0; j < dict_size(g->genes); ++j)
            count[g->covs[j][i]]++;
        for (j = 0; j < 101; ++j)
            fprintf(fp,"%s\t%d\t%d\n", dict_name(g->bcodes, i), j, count[j]);
    }

    fclose(fp);
    return 0;
}

int gene_cov(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    htsFile *fp  = hts_open(args.input_fname, "r");
    if (fp == NULL)
        error("%s : %s.", args.input_fname, strerror(errno));
    
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    hts_idx_t *idx = sam_index_load(fp, args.input_fname);
    if (idx == NULL)
        error("BAM file need be indexed.");

    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");

    hts_set_threads(fp, args.file_th);

    struct gtf_spec *G = gtf_read(args.gtf_fname, 1);
    if (G == NULL) error("GTF is empty.");

    struct gcov *gcov = gcov_init();
    // struct dict *CB_dict = dict_init();
    if (dict_read(gcov->bcodes, args.cell_list))
        error("Failed to load cell barcode.");
    
    struct dict *gene_dict = NULL;
    if (args.gene_list) {
        gene_dict = dict_init();
        if (dict_read(gene_dict, args.gene_list))
            error("Failed to load gene list.");
    }
        
    int i;
    bam1_t *b = bam_init1();
   
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *gl = &G->gtf[i];
        if (gene_dict) {
            int ret = dict_query(gene_dict, dict_name(G->gene_name, gl->gene_name));
            if (ret == -1) continue;
        }

        int tid;
        tid = bam_name2id(hdr, dict_name(G->name, gl->seqname));
        if (tid == -1) continue;
        int r;
        hts_itr_t *itr = sam_itr_queryi(idx, tid, gl->start, gl->end);
        // debug_print("Gene : %s", dict_name(G->gene_name, gl->gene_name));
        // each block come from exactly one gene        
        while ((r = sam_itr_next(fp, itr, b)) >= 0) {
            uint8_t *tag = bam_aux_get(b, args.tag);
            if (!tag) continue;
            char *cb = (char*)(tag+1);
            int ret;
            ret = dict_query(gcov->bcodes, cb);
            if (ret == -1) continue;
            push_gene_cov(gcov, dict_name(G->gene_name, gl->gene_name), cb, b);               
        }
        // count coverage of this block and reset buffer
        calc_buf_gene_cov(gcov, G, gl);
    }

    bam_destroy1(b);
    write_cov_mtx(args.output_fname, gcov);
    write_summary(args.summary_fname, gcov);

    gcov_destory(gcov);

    if (gene_dict) dict_destroy(gene_dict);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    gtf_destory(G);
    hts_close(fp);
    
    return 0;
}
