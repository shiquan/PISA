// count reads/fragments matrix for single-cell datasets
#include "utils.h"
#include "number.h"
#include "barcode_list.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

KHASH_MAP_INIT_STR(name, int)

typedef uint32_t count_t;

static int usage()
{
    fprintf(stderr, "* Count reads or fragments matrix for single-cell datasets.\n");
    fprintf(stderr, "Usage :\n  CountMatrix [options] aln.bam\n");
    fprintf(stderr, "Mandatory options :\n");
    fprintf(stderr, "    -tag        [CB]    Specify cell barcode tag.\n");
    fprintf(stderr, "    -anno_tag   [GN]    Annotation attribute.\n");
    fprintf(stderr, "    -list       [FILE]  Barcode list, white list, used as column names at matrix.\n");
    fprintf(stderr, "    -o          [FILE]  Output matrix.\n");
    fprintf(stderr, "    -umi        [UY]    UMI tag. Count once if more than one record has same UMI which overlapped with a region.\n");
    fprintf(stderr, "    -dis_corr            Disable correct UMI. Default all UMIs with 1 mismatch distance to each other are collapsed\n");
    fprintf(stderr, "    -q          [INT]   Minimal map quality to filter. [20]\n");
    fprintf(stderr, "    -count      [FILE]  UMI,Gene per cell barcode.\n");
    fprintf(stderr,"\n");

    return 1;
}
struct cell_barcode_counts {
    uint32_t nUMI;
    uint32_t nGene;
};

static struct args {
    const char *input_fname;
    const char *bed_fname;
    const char *tag; // cell barcode tag
    const char *anno_tag;
    const char *whitelist_fname;
    const char *output_fname;
    const char *umi_tag;
    const char *count_fname; // umi per cell barcode
    int mapq_thres;
    int dis_corr_umi;
    struct cell_barcode_counts *CBC;
} args = {
    .input_fname = NULL,
    .bed_fname = NULL,
    .tag = NULL,
    .anno_tag = NULL,
    .whitelist_fname = NULL,
    .output_fname = NULL,
    .umi_tag = NULL,
    .mapq_thres = 20,
    .dis_corr_umi = 0,
    .count_fname = NULL,
    .CBC = NULL,
};

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-anno_tag") == 0) var = &args.anno_tag;
        else if (strcmp(a, "-list") == 0) var = &args.whitelist_fname;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-dis_corr") == 0) {
            args.dis_corr_umi = 1;
            continue;
        }
        else if (strcmp(a, "-count") == 0) var = &args.count_fname;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument, %s", a);
    }
    if (args.input_fname == 0) error("No input bam.");
    if (args.output_fname == 0) error("No output file specified.");
    if (args.tag == 0) error("No cell barcode specified.");
    if (args.anno_tag == 0) error("No anno tag specified.");
    if (args.whitelist_fname == 0) error("No barcode list specified.");
    return 0;
}
struct mtx_counts {
    // int *v;
    int c; // counts
    int n, m;
    char **bcodes;
    kh_name_t *uhash;
};
static struct mtx_counts *mtx_counts_arr_init(int n)
{
    struct mtx_counts *m = malloc(n*sizeof(*m));
    int i;
    for (i = 0; i < n; ++i) {
        struct mtx_counts *m0 = &m[i];
        memset(m0, 0, sizeof(*m0));
        m0->uhash = kh_init(name);
    }
    return m;
}
static void mtx_counts_destory(struct mtx_counts *m, int n)
{    
    int i,j;
    for (i = 0; i < n; ++i) {
        struct mtx_counts *m0 = &m[i];
        for (j = 0; j < m0->n; ++j) free(m0->bcodes[j]);
        if (m0->bcodes) free(m0->bcodes);
        kh_destroy(name, m0->uhash);
    }
    free(m);    
}
static int check_similar(char *a, char *b)
{
    int l1, l2;
    l1 = strlen(a);
    l2 = strlen(b);
    if (l1 == 0 || l2 == 0) error("Try to compare an empty string");
    if (l1 != l2) error("Try to compare two unequal string.");
    int i, m = 0;
    for (i = 0; i < l1; ++i) {
        if (a[i] != b[i]) m++;
        if (m > 1) break;
    }
    if (m > 1) return 1;
    return 0;
}
static int get_val(struct mtx_counts *m, char *str)
{
    khint_t k = kh_get(name, m->uhash, str);
    assert (k != kh_end(m->uhash));
    return kh_val(m->uhash, k);
}
static void update_counts(struct mtx_counts **m, int l, int n)
{
    int i,j;
    for (i = 0; i < l; ++i) {
        struct mtx_counts *m0 = m[i];
        for (j = 0; j < n; ++j) {
            struct mtx_counts *m1 = &m0[j];
            if (args.dis_corr_umi == 1) {
                m1->c = m1->n;
                continue;
            }
            int *flag = malloc(m1->n*sizeof(int));
            memset(flag, 0, m1->n*sizeof(int));
            int i0, i1;
            for (i0 = 0; i0 < m1->n; ++i0) {
                if (flag[i0] == 1) continue;
                for (i1 = i0 + 1; i1 < m1->n; ++i1) {
                    if (flag[i1] == 1) continue;
                    if (check_similar(m1->bcodes[i0], m1->bcodes[i1]) == 0) {
                        int v0 = get_val(m1, m1->bcodes[i0]);
                        int v1 = get_val(m1, m1->bcodes[i1]);
                        if (v0 > v1) flag[i1] = 1;
                        else flag[i0] = 1;                        
                    }
                }
            }
            for (i0 = 0; i0 < m1->n; ++i0) {                
                if (flag[i0] == 0) m1->c++;
                // debug_print("%s\t%d", m1->bcodes[i0], flag[i0]);
            }
            free(flag);
        }
    }
}
int rank_cmp(const void *va, const void *vb)
{
    int a = *((const int*)va);
    int b = *((const int*)vb);
    return args.CBC[a].nUMI > args.CBC[b].nUMI;
              
}
int count_matrix(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return usage();

    struct lbarcode *lb = barcode_init();
    if (barcode_read(lb, args.whitelist_fname)) error("Empty white list.");

    htsFile *fp = hts_open(args.input_fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");
    
    FILE *out = fopen(args.output_fname, "w");
    CHECK_EMPTY(out, "%s : %s.", args.output_fname, strerror(errno));

    FILE *count = NULL;
    if (args.count_fname) {
        count = fopen(args.count_fname, "w");
        CHECK_EMPTY(count, "%s : %s.", args.count_fname, strerror(errno));
        args.CBC = malloc(lb->n*sizeof(struct cell_barcode_counts));
        memset(args.CBC, 0, sizeof(struct cell_barcode_counts)*lb->n);
    }
    
    kh_name_t *hash = kh_init(name); 
    int n=0, m=100;
    char **reg = malloc(m*sizeof(char*));
    struct mtx_counts **v = malloc(sizeof(void*)*m);
    
    bam1_t *b;
    bam1_core_t *c;
    int ret;
    b = bam_init1();
    c = &b->core;
    int id;
    khint_t k;
    
    while ((ret = sam_read1(fp, hdr, b)) >= 0) {
       
        uint8_t *tag = bam_aux_get(b, args.tag);
        if (!tag) error("Tag %s not found at line %s:%d\n", args.tag, hdr->target_name[c->tid], c->pos+1);
        id = barcode_select(lb, (char*)(tag+1));
        if (id == -1) continue;
        
        uint8_t *anno_tag = bam_aux_get(b, args.anno_tag);
        if (!anno_tag) continue;

        // Gene or Region
        char *val = (char*)(anno_tag+1);
        k = kh_get(name, hash, val);
        int r;
        if (k == kh_end(hash)) {
            if (m == n) {
                m = m *2;
                reg = realloc(reg, m*sizeof(char*));
                v = realloc(v, m*sizeof(int*));
            }
            reg[n] = strdup(val);
            k = kh_put(name, hash, reg[n], &r);
            kh_val(hash, k) = n;
            // v[n] = calloc(lb->n, sizeof(int));
            // memset(v[n], 0, lb->n*sizeof(int));
            v[n] = mtx_counts_arr_init(lb->n);
            n++;
        }
        
        int row = kh_val(hash, k);
        if (args.umi_tag) {
            uint8_t *umi_tag = bam_aux_get(b, args.umi_tag);
            if (!umi_tag) error("No UMI tag found at record. %s:%d", hdr->target_name[c->tid], c->pos+1);
            char *val = (char*)(umi_tag+1);
            struct mtx_counts *t = &v[row][id];
            k = kh_get(name, t->uhash, val);
            if (k == kh_end(t->uhash)) {
                if (t->n == t->m) {
                    t->m = t->m + 10;
                    t->bcodes = realloc(t->bcodes, t->m*sizeof(char*));
                }
                t->bcodes[t->n] = strdup(val);
                k = kh_put(name, t->uhash, t->bcodes[t->n], &r);
                kh_val(t->uhash, k) = 0;
                t->n++;
                // t->v[id]++;
            }
            else {
                kh_val(t->uhash, k)++;
            }
        }
        else 
            v[row][id].c++;
    }
    if (args.umi_tag) 
        update_counts(v, n, lb->n);

    if (count != NULL) {
        int i, j;
        for (i = 0; i < lb->n; ++i) {
            struct cell_barcode_counts *C = &args.CBC[i];
            for (j = 0; j < n; ++j) {
                if (v[i][j].c) {
                    C->nUMI += v[i][j].c;
                    C->nGene++;
                }
            }
        }
        int *rank;
        rank = malloc((lb->n+1)*sizeof(int));
        for (i = 0; i < lb->n+1; ++i) rank[i] = i+1;

        qsort(rank, lb->n, sizeof(int), rank_cmp);
        
        fprintf(count, "CELL_BARCODE\tnUMI\tnGene\trank\n");
        for (i = 0; i < lb->n; ++i) 
            fprintf(count, "%s\t%d\t%d\t%d\n", lb->b[i].s, args.CBC[i].nUMI, args.CBC[i].nGene, rank[i]);
        free(rank);
        free(args.CBC);
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    if (count) fclose(count);
    
    if (ret != -1) warnings("Truncated file?");   

    // header
    int i, j;
    fputs("ID", out);
    for (j = 0; j < lb->n; ++j) 
        fprintf(out, "\t%s",lb->b[j].s);
    fputc('\n', out);
    for (i = 0; i < n; ++i) {
        fprintf(out, "%s", reg[i]);
        for (j = 0; j < lb->n; ++j)  fprintf(out, "\t%d", v[i][j].c);
        fputc('\n', out);
    }
    fclose(out);

    kh_destroy(name,hash);
    for (i = 0; i < n; ++i) {
        free(reg[i]);
        mtx_counts_destory(v[i], lb->n);
    }
    free(reg); free(v);
    barcode_destory(lb);
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
