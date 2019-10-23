#include "utils.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "bam_pool.h"
#include "number.h"

struct impute_tags {
    char *name;
    char **tags;
};

KHASH_MAP_INIT_STR(name, int)

struct impute_index {
    int n, m;
    struct impute_tags *tags;
    kh_name_t *dict;
};
static struct impute_tags *impute_index_retrieve(struct impute_index *idx, char *s)
{
    khint_t k;
    k = kh_get(name, idx->dict, s);
    if (k == kh_end(idx->dict)) return NULL;
    return &idx->tags[kh_val(idx->dict, k)];
}

static struct args {
    const char *input_fname;
    const char *output_fname;
    int n_impute;
    char **impute_tags;
    int n_block;
    char **block_tags;

    htsFile *in;
    htsFile *out;
    bam_hdr_t *hdr;
    int qual_thres;
    int keep_all; // default will filter unannotated reads

    struct impute_index *index;

    int file_th;
    int n_thread;
    int chunk_size;
    uint64_t imputed_records;
    
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .n_impute = 0,
    .impute_tags = NULL,
    .n_block = 0,
    .block_tags = NULL,
    .in = NULL,
    .out = NULL,
    .hdr = NULL,
    .qual_thres = 20,
    .keep_all = 0,
    .index =NULL,
    .file_th = 1,
    .n_thread = 1,
    .chunk_size = 1000000,
    .imputed_records = 0,
};
static void impute_index_destory(struct impute_index *idx)
{
    int i;
    for (i = 0; i < idx->n; ++i) {
        int j;
        for (j = 0; j < args.n_impute; ++j)
            if (idx->tags[i].tags[j]) free(idx->tags[i].tags[j]);
        free(idx->tags[i].tags);
        free(idx->tags[i].name);
    }
    free(idx->tags);
    kh_destroy(name,idx->dict);
    free(idx);
}

static int parse_args(int argc, char **argv)
{
    int i;
    const char *block_tags = NULL;
    const char *impute_tags = NULL;
    const char *file_thread = NULL;
    const char *thread = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-impute") == 0) var = &impute_tags;
        else if (strcmp(a, "-block") == 0 ) var = &block_tags;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-k") == 0) {
            args.keep_all = 1;
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
    
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    CHECK_EMPTY(block_tags, "-block must be set.");
    CHECK_EMPTY(impute_tags, "-impute must be set.");
    
    kstring_t str = {0,0,0};
    kputs(block_tags, &str);
    int *s = ksplit(&str, ',', &args.n_block);
    assert(args.n_block>0);
    args.block_tags = malloc(sizeof(char*)*args.n_block);
    int j;
    for (j = 0; j < args.n_block; ++j) {
        args.block_tags[j] = strdup(str.s+s[j]);
        if (strlen(args.block_tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
    }
    free(s);
    str.l = 0;
    kputs(impute_tags, &str);
    int *s0 = ksplit(&str, ',', &args.n_impute);
    assert(args.n_impute > 0);
    args.impute_tags = malloc(sizeof(char*)*args.n_impute);
    for (j = 0; j < args.n_impute; ++j) {
        args.impute_tags[j] = strdup(str.s+s0[j]);
        if (strlen(args.impute_tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
    }
    free(s0);
    free(str.s);
        
    if (file_thread) args.file_th = str2int((char*)file_thread);
    assert(args.file_th>0);
    if (thread) args.n_thread = str2int((char*)thread);
    
    args.in  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.in, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    args.hdr = sam_hdr_read(args.in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    //int n_bed = 0;
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    if (args.file_th > 1) hts_set_threads(args.in, args.file_th);
    // hts_set_threads(args.out, args.file_th);
    
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
    return 0;
}
static void memory_release()
{
    impute_index_destory(args.index);
    bam_hdr_destroy(args.hdr);
    sam_close(args.in);
    sam_close(args.out);
    int i;
    for (i = 0; i < args.n_impute; ++i) free(args.impute_tags[i]);
    free(args.impute_tags);
    for (i = 0; i < args.n_block; ++i) free(args.block_tags[i]);
    free(args.block_tags);    
}

struct impute_tags_pool {
    int n, m;
    struct impute_tags *tags;
};

static struct impute_index *impute_index_build()
{
    struct impute_index *i = malloc(sizeof(*i));
    memset(i, 0, sizeof(*i));
    i->dict = kh_init(name);
    return i;
}
static void *build_idx_read(void *_p)
{
    struct impute_tags_pool *raw = malloc(sizeof(*raw));
    memset(raw, 0, sizeof(*raw));
    struct bam_pool *p = (struct bam_pool *)_p;
    
    kstring_t str = {0,0,0};
    int i;    
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = &p->bam[i];
        str.l = 0;
        if (raw->n == raw->m) {
            raw->m = raw->m == 0 ? 1024 : raw->m<<1;
            raw->tags = realloc(raw->tags, sizeof(struct impute_tags)*raw->m);
            raw->tags[raw->n].tags = malloc(args.n_impute*sizeof(void*)); // allocated
        }
        int j;
        for (j = 0; j < args.n_block; ++j) {
            uint8_t *tag = bam_aux_get(b, args.block_tags[j]);
            if (!tag) { str.l = 0; break; }
            kputs((char*)(tag+1), &str);
        }
        if (str.l == 0) continue;
        struct impute_tags *f = &raw->tags[raw->n];
        int is_empty = 1;        
        for (j = 0; j < args.n_impute; ++j) {
            uint8_t *tag = bam_aux_get(b, args.impute_tags[j]);
            f->tags[j] = NULL;
            if (!tag) continue;
            is_empty = 0;
            f->tags[j] = strdup((char*)(tag+1));
        }
        
        if (is_empty) continue; 
        f->name = strdup(str.s);
        raw->n++;
        if (raw->n < raw->m) 
            raw->tags[raw->n].tags = malloc(args.n_impute*sizeof(void*)); // init next record;
    }
    if (str.m) free(str.s);
    bam_pool_destory(p);
    return raw;
}
static void build_idx_core(struct impute_index *idx, struct impute_tags_pool *raw)
{
    int i;
    for (i = 0; i < raw->n; ++i) {
        struct impute_tags *r = &raw->tags[i];
        assert(r->name);
        khint_t k;
        int ret;
        char *s = strdup(r->name);
        k = kh_put(name, idx->dict, s, &ret);
        if (ret) {
            if (idx->n == idx->m) {
                idx->m = idx->m == 0 ? 1024 : idx->m<<1;
                idx->tags = realloc(idx->tags, idx->m*sizeof(struct impute_tags));
            }
            struct impute_tags *f = &idx->tags[idx->n];
            f->name = s;
            f->tags = malloc(sizeof(void*)*args.n_impute);
            int j;
            for (j = 0; j < args.n_impute; ++j)
                if (r->tags[j] == NULL) f->tags[j] = NULL;
                else f->tags[j] = strdup(r->tags[j]);
            kh_val(idx->dict, k) = idx->n++;            
        }
        else {
            int itr = kh_val(idx->dict, k);
            struct impute_tags *f = &idx->tags[itr];
            int j;
            for (j = 0; j < args.n_impute; ++j)
                if (f->tags[j] == NULL && r->tags[j]) f->tags[j] = strdup(r->tags[j]);
            free(s);
        }
        free(r->name);
        int j;
        for (j = 0; j < args.n_impute; ++j)
            if (r->tags[j]) free(r->tags[j]);
        free(r->tags);
    }
    if (raw->n < raw->m) free(raw->tags[raw->n].tags); // because we pre-allocate memory for next record ..
    free(raw->tags);
    free(raw);
}

static struct impute_index *build_index(const char *fn)
{
    htsFile *fp = hts_open(fn, "r");
    assert(fp);
    bam_hdr_t *hdr = sam_hdr_read(fp);
    hts_set_threads(fp, args.file_th);

    struct impute_index *index = impute_index_build();
    
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    for (;;) {
        struct bam_pool *b = bam_pool_create();
        bam_read_pool(b, fp, hdr, args.chunk_size);
        
        if (b == NULL) break;
        if (b->n == 0) { free(b->bam); free(b); break; }
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, build_idx_read, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct impute_tags_pool *d = (struct impute_tags_pool*)hts_tpool_result_data(r);
                build_idx_core(index, d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct impute_tags_pool *d = (struct impute_tags_pool*)hts_tpool_result_data(r);
        build_idx_core(index, d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    bam_hdr_destroy(hdr);
    sam_close(fp);

    return index;
}

struct p_data {
    struct bam_pool *bam;
    uint64_t c;
};

static void *run_it(void *data)
{
    struct bam_pool *p = (struct bam_pool*)data;
    struct p_data *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < p->n; ++i) {
        str.l = 0;
        bam1_t *b = &p->bam[i];
        int j;
        for (j = 0; j < args.n_block; ++j) {
            uint8_t *tag = bam_aux_get(b, args.block_tags[j]);
            if (!tag) { str.l = 0; break; }
            kputs((char*)(tag+1), &str);
        }
        if (str.l == 0) continue;
        struct impute_tags *f = impute_index_retrieve(args.index, str.s);
        if (f == NULL) continue;
        int set = 0;
        for (j = 0; j < args.n_impute; ++j) {
            uint8_t *tag = bam_aux_get(b, args.impute_tags[j]);
            if (!tag) {
                if (f->tags[j] == NULL) continue;
                // TODO: support type now
                bam_aux_append(b, args.impute_tags[j], 'Z', strlen(f->tags[j])+1, (uint8_t*)f->tags[j]);
                set = 1;
            }
        }
        if (set) d->c++;        
    }
    if (str.m) free(str.s);
    d->bam = p;
    return d;
}

static void write_out(struct p_data *p)
{
    args.imputed_records += p->c;
    int i;
    struct bam_pool *bam = p->bam;
    for (i = 0; i < bam->n; ++i)
        if (sam_write1(args.out, args.hdr, &bam->bam[i]) == -1) error("Failed to write SAM.");

    bam_pool_destory(bam);
    free(p);
}

// reads in each block with same BCs will be intepret as fragments from same template
// this program used to impute missed tags for reads from same block.
static int usage()
{
    fprintf(stderr, "LFR_impute in.bam\n");
    fprintf(stderr, "  -impute         Tags to impute.\n");
    fprintf(stderr, "  -block        Tags to identify each block.\n");
    fprintf(stderr, "  -k            Keep unclassified reads in the output.\n");
    fprintf(stderr, "  -@            Threads to pack and unpack bam file.\n");
    fprintf(stderr, "  -t            Threads to process.\n");
    return 1;
}

int LFR_impute(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();
    
    LOG_print("Build index ..");
    double t_real;
    t_real = realtime();        
    args.index = build_index(args.input_fname);
    
    LOG_print("Index build finished: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    for (;;) {
        struct bam_pool *b = bam_pool_create();
        bam_read_pool(b, args.in, args.hdr, args.chunk_size);
            
        if (b == NULL) break;
        if (b->n == 0) { free(b->bam); free(b); break; }
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct p_data *d = (struct p_data*)hts_tpool_result_data(r);
                write_out(d);   
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct p_data *d = (struct p_data*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    memory_release();    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    LOG_print("%"PRIu64" records updated.", args.imputed_records);
    return 0;
}
