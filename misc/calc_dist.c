#include "utils.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8193);

struct binRow {
    uint64_t *x;
    int cnt;
};
    
struct binMatrix {
    char **colnames; // column names
    char **rownames; // row names
    int n_col; // count of columns
    int n_row, m_row;
    int n; // array length
    struct binRow *r;
};

static void binary_matrix_destroy(struct binMatrix *M)
{
    int i;
    for (i = 0; i < M->n_row; ++i) free(M->rownames[i]);
    free(M->rownames);

    for (i = 0; i < M->n_col; ++i) if (M->colnames && M->colnames[i]) free(M->colnames[i]);
    if (M->colnames) free(M->colnames);

    for (i = 0; i < M->n_row; ++i) free(M->r[i].x);
    free(M->r);
    free(M);
}
static struct binMatrix *binary_matrix_load(const char *fname, int header, int transpose)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));

    struct binMatrix *M = malloc(sizeof(*M));
    memset(M, 0, sizeof(*M));
    
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int line = 0;
    int n;
    // int first_n = 0;
    int has_header = header;

    while (ks_getuntil(ks, 2, &str, &ret) >= 0) {

        line++;

        if (str.s[0] == '#') continue;
        if (str.l == 0) {
            warnings("Skip empty line. %d", line);
            continue;
        }

        int *s = ksplit(&str, '\t', &n);        
        if (M->n_col == 0 ) {
            M->n_col = n-1;
            M->n = M->n_col/64+1;
        }
        else if (M->n_col != n-1) error("Inconsistant columns at matrix. line %d.", line);
        
        if (has_header) {
            has_header = 0;
            M->colnames = malloc(sizeof(char*)*M->n_col);
            int i;
            for (i = 1; i < n; ++i) M->colnames[i-1] = strdup(str.s + s[i]);            
            continue;
        }

        if (M->n_row == M->m_row) {
            M->m_row = M->m_row == 0 ? 1024 : M->m_row*2;
            M->rownames = realloc(M->rownames, sizeof(char*)*M->m_row);
            M->r = realloc(M->r, sizeof(struct binRow)*M->m_row);
        }
        M->rownames[M->n_row] = strdup(str.s+s[0]);
        struct binRow *r = &M->r[M->n_row];
        memset(r, 0, sizeof(struct binRow));
        r->x = malloc(M->n*sizeof(uint64_t));
        memset(r->x, 0, M->n*sizeof(uint64_t));
        int i;
        for (i = 1; i < n; ++i) {
            int x = str2int(str.s+s[i]);
            if (x) {
                // Do not check the boundary, so it is may not be safe
                r->x[i/64] |= (1LL<<(i%64));
                r->cnt++;
            }
        }
        M->n_row++;
        free(s);
    }
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    if (M->n_row == 0 || M->n_col == 0) {
        free(M);
        return NULL;
    }
    return M;
}
static int usage()
{
    fprintf(stderr, "calc_dist input_matrix.txt\n");
    fprintf(stderr, " -t      [5]          Threads.\n");
    fprintf(stderr, " -method [Jaccard]    Distance method. Only support Jaccard now.\n");
    fprintf(stderr, " -o      [FILE]       Output matrix.\n");
    fprintf(stderr, " -r                   Transpose input matrix.\n");
    fprintf(stderr, " -h                   Skip first line. If -r set, transpose matrix first.\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *method;
    int n_thread;
    struct binMatrix *M;
    int transpose;
    int has_header;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .method = "Jaccard",
   
    .n_thread = 5,    
    .M = NULL,
    .transpose = 0,
    .has_header = 0,
};
static int parse_args(int argc, char **argv)
{
    if (argc == 1 ) return 1;
    const char *thread = NULL;
    int i;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-t") ==0) var = &thread;
        else if (strcmp(a, "-method") == 0) var = &args.method;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-r") == 0) {
            args.transpose = 1;
            continue;
        }
        else if (strcmp(a, "-h") == 0) {
            args.has_header = 1;
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

    CHECK_EMPTY(args.input_fname, "No input matrix.");
    CHECK_EMPTY(args.output_fname, "-o must be set.");

    if (thread) args.n_thread = str2int((char*)thread);
    if (args.n_thread < 1) args.n_thread = 1;
    
    args.M = binary_matrix_load(args.input_fname, args.has_header, args.transpose);
    CHECK_EMPTY(args.M, "Failed to load matrix");

    return 0;
}
static void write_matrix(float **V)
{
    int i, j;
    struct binMatrix *M = args.M;
    int n_row = M->n_row;

    FILE *fp_out = fopen(args.output_fname, "w");
    CHECK_EMPTY(fp_out,"%s : %s.", args.output_fname, strerror(errno));

    // header
    fputs("ID", fp_out);
    for (i = 0; i <n_row; ++i) {
        fputc('\t', fp_out);
        fputs(M->rownames[i], fp_out);
    }
    fputc('\n', fp_out);

    for (i = 0; i < n_row; i++) {        
        fputs(M->rownames[i], fp_out);

        for (j = 0; j < n_row; ++j) {
            if (i == j) V[i][j] = 0.5;
            if (i > j) V[i][j] = V[j][i];
            fputc('\t', fp_out);
            fprintf(fp_out, "%.2f", V[i][j]);
        }
        fputc('\n', fp_out);
    }

    fclose(fp_out);
}
static int countBits64(uint64_t x)
{
    int i;
    int c = 0;
    for (i = 0; i < 64; ++i) {
        if(x&(1LL<<i)) c++;
    }
    return c;
}

struct ret_v {
    int i;
    float *v;
};
static float calc_dist_jaccard(struct binMatrix *M, int i, int j)
{
    int k;
    struct binRow *x = &M->r[i];
    struct binRow *y = &M->r[j];
    if (x->cnt == 0 || y->cnt ==0) return 0.0;
    int ol = 0;
    for (k = 0; k < M->n; ++k) {
        if (x->x[k] > 0 && y->x[k] > 0) {
            uint64_t c = x->x[k] & y->x[k];
            if (c > 0) ol += countBits64(c);             
        }
    }
    return (float)ol/(x->cnt + y->cnt);
}
static void *run_it(void *_d)
{
    struct ret_v *d = (struct ret_v *)_d;
    struct binMatrix *M = args.M;
    float *v = malloc(M->n_row*sizeof(float));
    memset(v, 0, sizeof(float)*M->n_row);
    int i;    
    for (i = d->i+1; i < M->n_row; ++i) 
        v[i] = calc_dist_jaccard(M, i, d->i);
    d->v = v;
    return (void*)d;
}
int main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();
    

    struct binMatrix *M = args.M;
    
    float **V;
    V = malloc(M->n_row*sizeof(float*));
    memset(V, 0, sizeof(float*)*M->n_row);
    //V[0] = malloc(sizeof(float)*M->n_row);
    //memset(V[0], 0, sizeof(float)*M->n_row);
    
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    int i;
    for (i = 0; i < M->n_row; ++i) {
        int block;
        do {
            struct ret_v *d = malloc(sizeof(*d));
            d->i = i;
            block = hts_tpool_dispatch2(p, q, run_it, d, 1);
            if ((r = hts_tpool_next_result(q))) {                
                struct ret_v *d = (struct ret_v*)hts_tpool_result_data(r);
                V[d->i] = d->v;
                hts_tpool_delete_result(r, 1);
            }
        }
        while (block == -1);
    }
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct ret_v *d = (struct ret_v*)hts_tpool_result_data(r);
        V[d->i] = d->v;
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    write_matrix(V);
    binary_matrix_destroy(M);
    for (i = 0; i < M->n_row; ++i) free(V[i]);
    free(V);
    return 0;
}
