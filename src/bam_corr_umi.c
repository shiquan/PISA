#include "utils.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "umi_corr.h"
#include "bam_pool.h"
#include "string.h"
#include "htslib/kstring.h"
#include "number.h"
#include "htslib/thread_pool.h"

static int usage()
{
    fprintf(stderr, "bam_tag_corr in.bam\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o             Output bam.\n");
    fprintf(stderr, "  -tag           Tag for kmers to correct. Usually be UMI tag.\n");
    fprintf(stderr, "  -tags-block    Tags for each block. Reads in one block will be corrected with each other. Usually be cell barcode tag and gene tag.\n");
    fprintf(stderr, "  -@             Thread to unpack and pack BAM file.[5]\n");
    fprintf(stderr, "  -p             Thread to process data.[4]\n");
    // fprintf(stderr, "  -white-list    White list. Conflict with -block. If set tag will be correct with white list.\n");
    return 1;
}

static struct args {
    const char *input_fname;
    const char *output_fname;

    const char *tag;
    
    int n_block;
    char **blocks;
    
    htsFile *in;
    htsFile *out;

    bam_hdr_t *hdr;
    
    int n_thread;
    int file_th;
    struct corr_tag *Cindex;
    int update_count;

    int chunk_size;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .tag = NULL,
    .n_block = 0,
    .blocks = NULL,
    .hdr = NULL,
    .n_thread = 5,
    .file_th = 4,
    .Cindex = NULL,
    .update_count = 0,
    .chunk_size = 1000000, //1M
};

static int parse_args(int argc, char **argv)
{
    if (argc ==1) return 1;
    
    const char *block_tags = NULL;
    const char *file_th = NULL;
    const char *thread = NULL;
    
    int i;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-tags-block") == 0) var = &block_tags;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-p") == 0) var = &thread;

        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1] != '\0') error("Unknown option, %s", a);
        
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument : %s",a);
    }

    CHECK_EMPTY(args.output_fname, "-o need to be set.");
    CHECK_EMPTY(args.input_fname, "No input bam.");
    CHECK_EMPTY(args.tag, "-tag need to be set.");
    CHECK_EMPTY(block_tags, "-tags-block need to be set.");

    kstring_t str = {0,0,0};
    kputs(block_tags, &str);
    int *s = ksplit(&str, ',', &args.n_block);
    assert(args.n_block >0);
    args.blocks = malloc(args.n_block*sizeof(char*));
    for (i = 0; i < args.n_block; ++i)
        args.blocks[i] = strdup(str.s+s[i]);
    free(s);
    free(str.s);

    if (file_th) args.file_th = str2int((char*)file_th);
    if (thread) args.n_thread = str2int((char*)thread);

    args.in  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.in, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    args.hdr = sam_hdr_read(args.in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    hts_set_threads(args.in, args.file_th);
    hts_set_threads(args.out, args.file_th);
    
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
     
    return 0;
}

struct idx_str {
    char *key;
    char *val;
};
struct idx_str_pool {
    int n, m;
    struct idx_str *is;
};

static void build_idx_core(struct corr_tag *C, struct idx_str_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        struct idx_str *s = &p->is[i];
        corr_tag_push(C, s->key, s->val);
        free(s->key), free(s->val);
    }
    free(p->is);
    free(p);
}
static void *build_idx_read(void *data)
{
    struct bam_pool *d = (struct bam_pool*)data;
    struct idx_str_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < d->n; ++i) {
        if (p->n == p->m) {
            p->m = p->m == 0 ? 1024 : p->m<<1;
            p->is = realloc(p->is, sizeof(struct idx_str)*p->m);
        }
        bam1_t *b = &d->bam[i];
        str.l = 0;
        int j;
        for (j = 0; j < args.n_block; ++j) {
            uint8_t *tag = bam_aux_get(b, args.blocks[j]);
            if (!tag) { str.l = 0; break; }
        }
        if (str.l == 0) continue;
        uint8_t *tag = bam_aux_get(b, args.tag);
        if (!tag) continue;
        p->is[p->n].key = strdup(str.s);
        p->is[p->n].val = strdup((char*)(tag+1));
    }
    if (str.m) free(str.s);

    bam_pool_destory(d);
    return p;
}
static struct corr_tag *build_index(const char *fn)    
{
    htsFile *fp = hts_open(fn, "r");
    assert(fp);
    bam_hdr_t *hdr = sam_hdr_read(fp);
    hts_set_threads(fp, args.file_th);

    struct corr_tag *Cindex = corr_tag_build();
    
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
            block = hts_tpool_dispatch2(p, q, build_idx_read, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct idx_str_pool *d = (struct idx_str_pool*)hts_tpool_result_data(r);
                build_idx_core(Cindex, d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct idx_str_pool *d = (struct idx_str_pool*)hts_tpool_result_data(r);
        build_idx_core(Cindex, d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    bam_hdr_destroy(hdr);
    sam_close(fp);

    return Cindex;
}
struct p_data {
    struct bam_pool *p;
    int c;
};

static void *run_it(void *data)
{
    struct bam_pool *p = (struct bam_pool*)data;
    kstring_t str = {0,0,0};
    int i;
    int c = 0;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = &p->bam[i];
        str.l = 0;
        int j;
        for (j = 0; j < args.n_block; ++j) {
            uint8_t *tag = bam_aux_get(b, args.blocks[j]);
            if (!tag) { str.l = 0; break; }
        }
        if (str.l == 0) continue;
        uint8_t *tag = bam_aux_get(b, args.tag);
        if (!tag) continue;
        char *old_tag = (char*)(tag+1);
        char *new_tag = corr_tag_retrieve(args.Cindex, str.s, old_tag);
        if (strcmp(old_tag, new_tag) == 0) continue;
        if (bam_aux_update_str(b, args.tag, strlen(new_tag), new_tag))
            warnings("Failed to update tag. %s", (char*)b->data);
        c++;
    }

    struct p_data *pd = malloc(sizeof(*pd));
    pd->p = p;
    pd->c = c;
    return pd;
}
static void write_out(struct p_data *pd)
{
    struct bam_pool *p = pd->p;
    int i;
    for (i = 0; i < p->n; ++i)        
        if (sam_write1(args.out, args.hdr, &p->bam[i]) == -1) error("Failed to write SAM.");
    args.update_count += pd->c;
    bam_pool_destory(p);
    free(pd);
}
static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.in);
    sam_close(args.out);
    int i;
    for (i = 0; i < args.n_block; ++i) free(args.blocks[i]);
    free(args.blocks);
    corr_tag_destory(args.Cindex);
}

int bam_corr_umi(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    LOG_print("Build index ..");
    double t_real;
    t_real = realtime();        
    args.Cindex = build_index(args.input_fname);
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
    LOG_print("%d records updated.", args.update_count);
    return 0;
}
