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
    fprintf(stderr, "    -tag        Specify cell barcode tag.\n");
    fprintf(stderr, "    -anno_tag   Annotation attribute.\n");
    fprintf(stderr, "    -list       Barcode list, white list, used as column names at matrix.\n");
    fprintf(stderr, "    -o          Output matrix.\n");
    fprintf(stderr, "    -umi        UMI tag. Count once if more than one record has same UMI which overlapped with a region.\n");
    // fprintf(stderr, "    -t          Threads. [5]\n");
    fprintf(stderr, "    -q          Minimal map quality to filter. [60]\n");
    // fprintf(stderr, "    -i          Maximal insert size to filter, 0 for not set. [0]\n");
    fprintf(stderr,"\n");

    return 1;
}

static struct args {
    const char *input_fname;
    const char *bed_fname;
    const char *tag; // cell barcode tag
    const char *anno_tag;
    const char *whitelist_fname;
    const char *output_fname;
    const char *umi_tag;
    int mapq_thres;
    
} args = {
    .input_fname = NULL,
    .bed_fname = NULL,
    .tag = NULL,
    .anno_tag = NULL,
    .whitelist_fname = NULL,
    .output_fname = NULL,
    .umi_tag = NULL,
    .mapq_thres = 60,
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
    int *v;
    int n, m;
    char **bcodes;
    kh_name_t *uhash;
};
static struct mtx_counts *mtx_counts_init(int n)
{
    struct mtx_counts *m = malloc(sizeof(*m));
    memset(m, 0, sizeof(*m));
    m->v = calloc(n, sizeof(int));
    m->uhash = kh_init(name);    
    return m;
}
static void mtx_counts_destory(struct mtx_counts *m)
{
    int i;
    for (i = 0; i < m->n; ++i) free(m->bcodes[i]);
    if (m->bcodes) free(m->bcodes);
    kh_destroy(name, m->uhash);
    free(m->v);
    free(m);    
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
            v[n] = mtx_counts_init(lb->n);

            n++;
        }
        int row = kh_val(hash, k);
        if (args.umi_tag) {
            uint8_t *umi_tag = bam_aux_get(b, args.umi_tag);
            if (!umi_tag) error("No UMI tag found at record. %s:%d", hdr->target_name[c->tid], c->pos+1);
            char *val = (char*)(umi_tag+1);
            struct mtx_counts *t = v[row];
            k = kh_get(name, t->uhash, val);
            if (k == kh_end(t->uhash)) {
                if (t->n == t->m) {
                    t->m = t->m + 10;
                    t->bcodes = realloc(t->bcodes, t->m*sizeof(char*));
                }
                t->bcodes[t->n] = strdup(val);
                k = kh_put(name, t->uhash, t->bcodes[t->n], &r);
                t->n++;
                t->v[id]++;
            }
        }
        else 
            v[row]->v[id]++;       
    }    
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    if (ret != -1) warnings("Truncated file?");
    int i, j;
    fputs("ID", out);
    for (j = 0; j < lb->n; ++j) 
        fprintf(out, "\t%s",lb->b[j].s);
    fputc('\n', out);
    for (i = 0; i < n; ++i) {
        fprintf(out, "%s", reg[i]);
        for (j = 0; j < lb->n; ++j)  fprintf(out, "\t%d", v[i]->v[j]);
        fputc('\n', out);
    }
    fclose(out);

    kh_destroy(name,hash);
    for (i = 0; i < n; ++i) {
        free(reg[i]);
        mtx_counts_destory(v[i]);
    }
    free(reg); free(v);
    barcode_destory(lb);
    return 0;
}
