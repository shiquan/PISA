#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "fastq.h"
#include "barcode_list.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include <zlib.h>
#include <ctype.h>
#include <sys/stat.h>

KSEQ_INIT(gzFile, gzread)
static int usage()
{
    fprintf(stderr, "* Sort fastq records by tags and deduplicate.\n");
    fprintf(stderr, "fastq-sort  in.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, " -tag       Tags, such as CB,UR. Order of these tags is sensitive.\n");
    fprintf(stderr, " -dedup     Remove dna copies with same tags. Only keep reads have the best quality.\n");
    fprintf(stderr, " -list      White list for first tag, usually for cell barcodes.\n");
    fprintf(stderr, " -t         Threads to deduplicate.\n");
    fprintf(stderr, " -o         Output fastq.");
    fprintf(stderr, " -p         Input fastq is smart pairing.\n");
    fprintf(stderr, " -dropN     Drop if N found in tags.\n");
    fprintf(stderr, " -mem       In memory mode. Put all records in memory.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options for splited output files:\n");
    fprintf(stderr, " -split     Split records to files by first tag. -O and -list must be set with this parameter.\n");
    fprintf(stderr, " -O         Outdir of split fastq. Conflict with -o.\n");
    fprintf(stderr, "Notice:\n");
    fprintf(stderr, "     1. This program is designed for single cell data. Deduplicate step is computing \n");
    fprintf(stderr, "        extensive if too much reads in one block. It is not suggest apply to bulk data.\n");
    fprintf(stderr, "     2. For index mode (default), input should be unpacked.\n");
    fprintf(stderr, "     3. In case you don't care about the memory and need the result ASAP. -mem may be a choice.\n");
    return 1;
}

struct data_block {
    long start;
    long end;
};

struct data_index;
KHASH_MAP_INIT_STR(idx, struct data_index*)

struct data_index {
    char *fkey; // point to key, don't free it
    int n, m;
    char **keys;
    int n_idx, m_idx;
    struct data_block *idx;
    struct bseq_pool *pool;
    kh_idx_t *dict;
};
struct data_index *data_index_create()
{
    struct data_index *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    return d;
}
static struct data_index *data_idx_retrieve(struct data_index *idx, char *s)
{
    khint_t k;
    if (idx->dict == 0) {
        idx->dict = kh_init(idx);
        goto idx_hash_push;
    }

    k = kh_get(idx, idx->dict, s);
    if (k != kh_end(idx->dict)) return kh_val(idx->dict,k);

  idx_hash_push:
    if (idx->n == idx->m) {
        idx->m = idx->m == 0 ? 10 : idx->m<<1;
        idx->keys = realloc(idx->keys, idx->m*sizeof(char*));
    }
    
    idx->keys[idx->n] = strdup(s);
    int ret;
    k = kh_put(idx, idx->dict, idx->keys[idx->n], &ret);
    struct data_index *i0 = data_index_create();
    i0->fkey = idx->keys[idx->n];
    kh_val(idx->dict, k) = i0;
    idx->n++;
    
    return kh_val(idx->dict, k);
}

static void data_index_clean(struct data_index *d)
{
    if (d->n > 0) {
        int i;
        for (i = 0; i < d->n; i++) {
            struct data_index *i0 = data_idx_retrieve(d, d->keys[i]);
            data_index_clean(i0);
            free(d->keys[i]);
            free(i0);
        }
        free(d->keys);
        kh_destroy(idx, d->dict);
    }
    else {
        if (d->m_idx) free(d->idx);
        if (d->pool) bseq_pool_destroy(d->pool);
    }
    memset(d, 0, sizeof(*d));
}
void data_index_destory(struct data_index *d)
{
    data_index_clean(d);
    free(d);
}

/*
struct FILE_tag_index {
    int n, m;
    char **key;
    void *dict;
};
*/
static struct args {
    const char *input_fname;
    const char *list_fname;
    const char *output_fname;
    const char *outdir;
    int split_file;
    int dropN;    
    int n_thread;
    int dedup;
    int n_tag;
    int in_mem;
    int smart_pairing;
    int check_list;

    int file_count;
    
    char **tags;
    FILE *out;
    // For each thread, keep file handler
    //FILE **fp;
    FILE *fp_in;
    struct lbarcode *lb;
    long end;
    
} args = {
    .input_fname = NULL,
    .list_fname = NULL,
    .output_fname = NULL,
    .outdir = NULL,
    .split_file = 0,
    .dropN = 0,
    .smart_pairing = 0,
    .check_list = 0,
    .file_count = 0,
    .dedup = 0,
    .n_tag = 0,
    .in_mem = 0,
    .tags = NULL,
    .out = NULL,
    .fp_in = NULL,
    .lb = NULL,
    .end = 0,
};

#define FILE_PER_FOLD 100
#define MAX_FILE_NUM  10000
static void write_out_bseq(struct bseq_pool *pool, FILE *fp)
{
    int i;
    for (i = 0; i < pool->n; ++i) {
        struct bseq *b = &pool->s[i];
        if (b->flag == FQ_PASS) {            
            fprintf(fp, "%s\n%s\n+\n%s\n", b->n0, b->s0, b->q0);
            if (b->l1) fprintf(fp, "%s\n%s\n+\n%s\n", b->n0, b->s1, b->q1);
        }
    }
    // bseq_pool_destroy(pool);
}
static void write_out_core(struct data_index *idx, FILE *fp)
{
    if (idx->n > 0) {
        int i;
        for (i = 0; i < idx->n; ++i) {
            struct data_index *i0 = data_idx_retrieve(idx, idx->keys[i]);
            write_out_bseq(i0->pool, fp);
        }
    }
    else
        write_out_bseq(idx->pool, fp);
}
void write_out(struct data_index *idx, int c)
{
    if (args.split_file) {
        int sub0 = c/FILE_PER_FOLD;
        kstring_t str = {0,0,0};
        if (args.outdir) kputs(args.outdir, &str); // / already added
        if (str.s[str.l-1] != '/') kputc('/', &str);
        kputw(sub0,&str);kputc('/', &str);
        mkdir(str.s, 0777);
        kputs(idx->fkey, &str);
        kputs(".fq", &str);
        FILE *fp = fopen(str.s, "w");
        CHECK_EMPTY(fp, "%s : %s.", str.s, strerror(errno));
        write_out_core(idx, fp);
        fclose(fp);
        free(str.s);
    }
    else {
        assert(args.out);
        write_out_core(idx, args.out);
    }
    data_index_clean(idx);
}
/*
static void FILE_tag_destroy(struct FILE_tag_index *idx)
{
    int i;
    for (i = 0; i < idx->n; ++i) {
        khiter_t k;
        k = kh_get(idx, (kh_idx_t*)idx->dict, idx->key[i]);
        struct data_index *di = kh_val((kh_idx_t*)idx->dict, k);
        free(di->idx);
        free(di);
        kh_del(idx,(kh_idx_t*)idx->dict,k);
        free(idx->key[i]);
    }
    free(idx->key);
    kh_destroy(idx,(kh_idx_t*)idx->dict);
    free(idx);
}
*/
struct bseq *bend_to_bseq(char *s)
{
    kstring_t str = {0,0,0};
    kputs(s, &str);
    int n;
    int *p = ksplit(&str, '\n', &n);
    if (n != 4 && n!= 8) {
        fprintf(stderr, "%d\n%s\n", n, s);
        assert(1);
    }
    struct bseq *b = malloc(sizeof(*b));
    memset(b, 0, sizeof(struct bseq));
    b->n0 = strdup(str.s+p[0]);
    b->s0 = strdup(str.s+p[1]);
    b->q0 = strdup(str.s+p[3]);
    b->l0 = strlen(b->s0);
    int l = strlen(b->q0);
    if (l != b->l0) error("Unequal sequence and quality length. %s vs %s", b->s0, b->q0);
    if (args.smart_pairing) {
        b->s1 = strdup(str.s+p[5]);
        b->q1 = strdup(str.s+p[7]);
        b->l1 = strlen(b->s1);
        l = strlen(b->q1);
        if (l != b->l1) error("Unequal sequence and quality length. %s vs %s", b->s1, b->q1);
    }
    
    free(str.s);
    return b;
}
static struct bseq *read_file_block(FILE *fp, struct data_block *block)
{
    int l = block->end - block->start;    
    if (fseek(fp, block->start,SEEK_SET)) error("failed to seek.");
    char *s = malloc(sizeof(char)*(block->end -block->start+1));
    if (s == NULL) error("Fail to allocate memory.");    
    fread(s, sizeof(char), block->end-block->start, fp);
    s[l] = '\0';
    kstring_t str = {0,0,0};
    kputsn(s, 100, &str);
    // debug_print("%llu\t%s", block->start, str.s);
    free(str.s);
    struct bseq *b = bend_to_bseq(s);
    free(s);
    return b;
}
static void test_read_file_block(FILE *fp, struct data_block *block)
{

    int l = block->end - block->start;
    
    fseek(fp, block->start,SEEK_SET);
    char *s = malloc(sizeof(char)*(block->end -block->start+1));    
    fread(s, sizeof(char), block->end-block->start, fp);
    s[l] = '\0';

    puts(s);

    puts("\n=======\n");
    
    free(s);
    // return NULL;
}
/*
void test_file_index(struct FILE_tag_index *idx, const char *fn)
{
    FILE *fp = fopen(fn, "r");
    
    int i, j;
    for (i = 0; i < idx->n; ++i) {
        khint_t k;
        k = kh_get(idx, (kh_idx_t*)idx->dict, idx->key[i]);
        struct data_index *di = kh_val((kh_idx_t*)idx->dict, k);
        assert(di);
        for (j = 0; j < di->n; ++j) {
            test_read_file_block(fp, &di->idx[j]);
        }
    }

    fclose(fp);
}

char *build_string_from_array(int n, char **v)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < n; ++i) kputs(v[i], &str);
    return str.s;
}
*/
static void memory_release()
{
    // fclose(args.out);
    int i;
    //for (i = 0; i < args.n_thread; ++i) fclose(args.fp[i]);
    fclose(args.fp_in);
    for (i = 0; i < args.n_tag; ++i) free(args.tags[i]);
    free(args.tags);
    if (args.check_list) barcode_destory(args.lb);
}
static int FILE_read_line_length(FILE *fp)
{
    int l;
    for (l= 0; fgetc(fp) != '\n';++l);
    return l;
}
static void push_idx_idx_core(struct data_index *idx, long start, long end)
{
    if (idx->n_idx == idx->m_idx) {
        idx->m_idx = idx->m_idx == 0 ? 10 : idx->m_idx<<1;
        idx->idx = realloc(idx->idx, idx->m_idx *sizeof(struct data_block));
    }
    idx->idx[idx->n_idx].start = start;
    idx->idx[idx->n_idx].end = end;
    idx->n_idx++;
}
void push_idx_idx(struct data_index *idx, int n, char **names, long start, long end)
{

    struct data_index *i0 = data_idx_retrieve(idx, names[0]);

    if (n > 1) {
        kstring_t str= {0,0,0};
        int i;
        for (i = 1; i < n;++i) kputs(names[i], &str);
        struct data_index *i1 = data_idx_retrieve(i0, str.s);
        push_idx_idx_core(i1, start, end);

        free(str.s);
    }
    else
        push_idx_idx_core(i0, start, end);
}
static void push_idx_mem_core(struct data_index *idx, struct bseq *b)
{
    if (idx->pool == NULL)
        idx->pool = bseq_pool_init();
    bseq_pool_push(b, idx->pool);
}
static void push_idx_mem(struct data_index *idx, int n, char **names, struct bseq *b)
{
    struct data_index *i0 = data_idx_retrieve(idx, names[0]);
    if (n > 1) {
        kstring_t str= {0,0,0};
        int i;
        for (i = 1; i < n;++i) kputs(names[i], &str);
        struct data_index *i1 = data_idx_retrieve(i0, str.s);
        push_idx_mem_core(i1, b);
        free(str.s);
    }
    else
        push_idx_mem_core(i0, b);
}
//static struct FILE_tag_index *build_file_index(const char *fname, int in_mem)
static struct data_index *build_file_index(const char *fname, int in_mem)
{
    struct data_index *idx = malloc(sizeof(*idx));
    memset(idx, 0, sizeof(*idx));
    // idx->dict = kh_init(idx);
        
    if (args.in_mem) {
        gzFile fp = gzopen(fname, "r");
        kseq_t *ks;
        ks = kseq_init(fp);
        struct bseq b = {0,0,0,0,0,0,0,0};
        while (kseq_read(ks)>=0) {
            memset(&b, 0, sizeof(b));
            b.n0 = strdup(ks->name.s);
            b.s0 = strdup(ks->seq.s);
            b.l0 = ks->seq.l;
            b.q0 = strdup(ks->qual.s);

            if (args.smart_pairing) {
                if (kseq_read(ks) < 0) error("Failed to found read 2.");
                if (strcmp(ks->name.s, b.n0) != 0) error("Inconsistant read name, %s vs %s", ks->name.s, b.n0);
                b.s1 = strdup(ks->seq.s);
                b.q1 = strdup(ks->qual.s);
                b.l1 = ks->seq.l;
            }
            char **names = fastq_name_pick_tags(ks->name.s, args.n_tag, args.tags);
            int i;
            if (args.check_list) {
                int ret = barcode_select(args.lb, names[0]);
                if (ret == -1) {
                    
                    for (i = 0; i < args.n_tag; ++i) free(names[i]);
                    free(names);
                    continue;  
                }
            }
                        
            push_idx_mem(idx,args.n_tag, names, &b);
            for (i = 0; i < args.n_tag; ++i) free(names[i]);
            free(names);
        }
        kseq_destroy(ks);
        gzclose(fp);
    }
    else {
        FILE *fp = fopen(fname, "r");
        CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));
        fseek(fp,0L,SEEK_END);
        fflush(fp);
        args.end = ftell(fp);
        // debug_print("%d", args.end);
        fseek(fp,0L, SEEK_SET);
        
        LOG_print("The file is %ld B.", args.end);
        // struct FILE_tag_index *idx = malloc(sizeof(*idx));

        kstring_t str = {0,0,0};
        kstring_t str2 = {0,0,0}; // for read 2, used in smart pairing mode
        long start, end;
        int begin_of_record = 1;
        
        for (;;) {
            int c = fgetc(fp);
            if (c == EOF) break;
            
            str.l = 0; // clear buffer
            
            if (begin_of_record) {
                fflush(fp);
                start = ftell(fp);
                begin_of_record = 0;
                assert (c == '@');
            }
            for (; !isspace(c) && c != '\n';) {
                kputc(c, &str);
                c = fgetc(fp);
            }
            if (isspace(c) && c!= '\n') {
                for (;c!= '\n' && c!= EOF; c= fgetc(fp)); // emit tails
                // c = fgetc(fp); // emit '\n'
            }
            
            char **names = fastq_name_pick_tags(str.s, args.n_tag, args.tags);  
            int l = FILE_read_line_length(fp);
            fseek(fp, l+3, SEEK_CUR);
            begin_of_record = 1;
            
            if (args.smart_pairing) { // the next read is paired end
                str2.l = 0; // clean buffer
                c = fgetc(fp);
                while (isspace(c)) c = fgetc(fp); // skip empty line
                assert(c == '@'); // read names
                for (;c != '\n';) {
                    c = fgetc(fp);
                    kputc(c, &str2);
                }
                if (strcmp(str.s, str2.s) != 0) error("Inconstant read name at smart pairing mode. %s vs %s.", str.s, str2.s);
                int l = FILE_read_line_length(fp);
                fseek(fp, l+3, SEEK_CUR);            
            }
            fflush(fp);
            end = ftell(fp);
            
            assert(end <= args.end);
            // block FILE:start-end store this read
            int i;
            if (names == NULL) continue; // no tags in the read names
            if (args.check_list) {
                //debug_print("%s", names[0]);
                int ret = barcode_select(args.lb, names[0]);
                if (ret == -1) {
                    for (i = 0; i < args.n_tag; ++i) free(names[i]);
                    free(names);
                    continue;  
                }
            }
            
            push_idx_idx(idx, args.n_tag, names, start, end);
        }
        fclose(fp);
        
    }
    return idx;
}
/*
static void *run_it(void *_d, int i)
{
    struct data_index *idx = (struct data_index*)_d;    
    struct bseq_pool *p = bseq_pool_init();

    // FILE *fp = args.fp[i];
    
    int k;   
    for (k = 0; k < idx->n; ++k) {
        struct bseq *b = read_file_block(fp, &idx->idx[k]);
        bseq_pool_push(b, p);
    }
    if (args.dedup == 1)
        bseq_pool_dedup(p);
    return p;        
}
*/
static void *run_it1(void *_d)
{
    /*
    struct bseq_pool *p = (struct bseq_pool*)_d;
    if (p->n > 1) 
        bseq_pool_dedup(p);
    */
    struct data_index *p = (struct data_index*)_d;
    if (args.dedup == 0) return p; // do nothing, here only dedup 
    // per tag
    
    if (p->n > 0) { // have keys, more than one tags            
        int i;
        for (i = 0; i < p->n; ++i) {    
            struct data_index *i0 = data_idx_retrieve(p, p->keys[i]);
            i0->fkey = p->keys[i];
            bseq_pool_dedup(i0->pool);
        }
    }
    else {
        bseq_pool_dedup(p->pool);
    }
    return p;        
}

//static void write_out(struct bseq_pool *p)
/*
static void write_out(struct data_index *p)
{

    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        if (b->flag == FQ_PASS) {            
            fprintf(args.out, "%s\n%s\n+\n%s\n", b->n0, b->s0, b->q0);
            if (b->l1) fprintf(args.out, "%s\n%s\n+\n%s\n", b->n0, b->s1, b->q1);
        }
    }
    bseq_pool_destroy(p);

}
*/
static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *thread = NULL;
    const char *tags = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-list") == 0) var = &args.list_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-O") == 0) var = &args.outdir;
        else if (strcmp(a, "-dropN") == 0) {
            args.dropN = 1;
            continue;
        }
        else if (strcmp(a, "-dedup") == 0) {
            args.dedup = 1;
            continue;
        }
        else if (strcmp(a, "-p") == 0) {
            args.smart_pairing = 1;
            continue;
        }
        else if (strcmp(a, "-mem") == 0) {
            args.in_mem = 1;
            continue;
        }
        else if (strcmp(a, "-split") == 0) {
            args.split_file = 1;
            continue;
        }

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
    // todo: check file, NO streaming allowed
    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (tags == NULL) error("-tag must be set.");

    kstring_t str = {0,0,0};
    kputs(tags, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    args.n_tag = n;
    args.tags = malloc(n*sizeof(char*));
    for (i = 0; i <n; ++i) args.tags[i] = strdup(str.s+s[i]);
    free(s); free(str.s);

    if (thread) args.n_thread = str2int((char*)thread);
    if (args.n_thread < 1) args.n_thread = 1;

    if (args.list_fname) {
        args.lb = barcode_init();
        if (barcode_read(args.lb, args.list_fname)) error("Barcode list is empty.");
        args.check_list = 1;
        if (args.split_file && args.lb->n > MAX_FILE_NUM)
            error("Too much split files. Try to reduce you list under %d", MAX_FILE_NUM);
    }

    if (args.split_file) {
        if (args.outdir == NULL) error("Please set -O when use -split.");
        if (args.output_fname) warnings("-o will be ignore when set -split.");
        if (args.check_list == 0) error("-list must be set if use -split.");
    }

    if (args.outdir == NULL && args.split_file == 0) {
        args.out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
        CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    }
    return 0;
}

int fill_pool_by_index(struct data_index *idx)
{
    if (args.in_mem == 1) return 0;
    if (idx->n > 0) {
        int i;
        for (i = 0; i < idx->n; ++i) {
            struct data_index *i0 = data_idx_retrieve(idx, idx->keys[i]);
            fill_pool_by_index(i0);
        }
    }
    else {
        assert(idx->pool == NULL);
        idx->pool = bseq_pool_init();
        int j;
        for (j = 0; j < idx->n_idx; ++j) {
            struct bseq *b = read_file_block(args.fp_in, &idx->idx[j]);
            bseq_pool_push(b, idx->pool);
            free(b);
        }        
    }
    return 0;
}
int fsort(int argc, char ** argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();

    //struct FILE_tag_index *idx = build_file_index(args.input_fname);
    struct data_index *idx = build_file_index(args.input_fname, args.in_mem);

    LOG_print("Build index time: %.3f sec; CPU: %.3f sec", realtime()-t_real, cputime());

    // thread cache
    /*
    args.fp = malloc(args.n_thread*sizeof(void*));
    int i;
    for (i = 0; i < args.n_thread; ++i) args.fp[i] = fopen(args.input_fname, "r");
    */

    // test_file_index(idx, args.input_fname);

    int i;
    // khint_t k;
    if (args.in_mem == 0) 
        args.fp_in = fopen(args.input_fname,"r");
    
    if (args.dedup) {
        hts_tpool *pool = hts_tpool_init(args.n_thread);
        hts_tpool_process *q = hts_tpool_process_init(pool, args.n_thread*2, 0);
        hts_tpool_result *r;        
        for (i = 0; i < idx->n; ++i) {
            struct data_index *i0 = data_idx_retrieve(idx, idx->keys[i]);
            fill_pool_by_index(i0);

            int block;
            do {
                block = hts_tpool_dispatch2(pool, q, run_it1, i0, 1);
                if ((r = hts_tpool_next_result(q))) {                
                    struct data_index  *di = (struct data_index*)hts_tpool_result_data(r);
                    write_out(di, args.file_count++);
                    hts_tpool_delete_result(r, 0);
                }
            }
            while (block == -1);
        }
        hts_tpool_process_flush(q);
        
        while ((r = hts_tpool_next_result(q))) {
            struct data_index *di = (struct data_index*)hts_tpool_result_data(r);
            write_out(di, args.file_count++);
            hts_tpool_delete_result(r, 0);
        }
        hts_tpool_process_destroy(q);
        hts_tpool_destroy(pool);
    }
    else {
    
        for (i = 0; i < idx->n; ++i) {
            struct data_index *i0 = data_idx_retrieve(idx, idx->keys[i]);
            fill_pool_by_index(i0);
            write_out(i0, args.file_count++);
        }
    }

    //FILE_tag_destroy(idx);
    data_index_destory(idx);
    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());

    return 0;
}

/*
int main(int argc, char **argv)
{
    return fsort(argc, argv);
}
*/
