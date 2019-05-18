#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "fastq.h"
#include "barcode_list.h"
#include "number.h"
#include "thread_pool.h"
#include <ctype.h>

static int usage()
{
    fprintf(stderr, "* Sort fastq records by tags and deduplicate. Input should be unpacked.\n");
    fprintf(stderr, "fastq-sort  in.fq\n");
    fprintf(stderr, " -tag       Tags, such as CB,UR. Order of these tags is sensitive.\n");
    fprintf(stderr, " -dedup     Remove dna copies with same tags. Keep reads with best quality.\n");
    fprintf(stderr, " -list      White list for first tag.\n");
    fprintf(stderr, " -t         Threads.\n");
    fprintf(stderr, " -o         Output fastq.");
    fprintf(stderr, " -p         Input fastq is smart pairing.\n");
    fprintf(stderr, " -dropN     Drop if N in tags.\n");
    fprintf(stderr, "Notice: This program is designed for single cell data. Deduplicate step is computing \n");
    fprintf(stderr, "        extensive if too much reads in one block. It is not suggest apply to bulk data.\n");
    return 1;
}

struct data_block {
    uint64_t start;
    uint64_t end;
};

struct data_index {
    int n, m;
    struct data_block *idx;
};

KHASH_MAP_INIT_STR(idx, struct data_index*)

struct FILE_tag_index {
    int n, m;
    char **key;
    void *dict;
};
static struct args {
    const char *input_fname;
    const char *list_fname;
    const char *output_fname;
    int dropN;    
    int n_thread;
    int dedup;
    int n_tag;
    int smart_pairing;
    int check_list;
    char **tags;
    FILE *out;
    // For each thread, keep file handler
    FILE **fp;
    struct lbarcode *lb;
    uint64_t end;
} args = {
    .input_fname = NULL,
    .list_fname = NULL,
    .output_fname = NULL,
    .dropN = 0,
    .smart_pairing = 0,
    .check_list = 0,
    .dedup = 0,
    .n_tag = 0,  
    .tags = NULL,
    .out = NULL,
    .fp = NULL,
    .lb = NULL,
    .end = 0,
};

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
static void memory_release()
{
    fclose(args.out);
    int i;
    for (i = 0; i < args.n_thread; ++i) fclose(args.fp[i]);
    free(args.fp);
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
static struct FILE_tag_index *build_file_index(const char *fname)
{
    FILE *fp = fopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));
    fseek(fp,0,SEEK_END);
    args.end = ftell(fp);
    // debug_print("%d", args.end);
    fseek(fp,0, SEEK_SET);
    struct FILE_tag_index *idx = malloc(sizeof(*idx));
    memset(idx, 0, sizeof(*idx));
    idx->dict = kh_init(idx);
    kstring_t str = {0,0,0};
    kstring_t str2 = {0,0,0}; // for read 2, used in smart pairing mode
    uint64_t start, end;
    int begin_of_record = 1;

    for (;;) {
        int c = fgetc(fp);
        if (c == EOF) break;
        
        str.l = 0; // clear buffer
        
        if (begin_of_record) {
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
        // push to index
        char *s = build_string_from_array(args.n_tag, names);
        for (i = 0; i < args.n_tag; ++i) free(names[i]);
        free(names);        
        khint_t k = kh_get(idx, idx->dict, s);
        if (k == kh_end((kh_idx_t*)idx->dict)) {
            if (idx->n == idx->m) {
                idx->m = idx->m == 0 ? 1024 : idx->m<<1;
                idx->key = realloc(idx->key, idx->m*sizeof(char*));
            }
            idx->key[idx->n]= strdup(s);
            int ret;
            k = kh_put(idx, idx->dict, idx->key[idx->n], &ret);
            struct data_index *di = malloc(sizeof(*di));
            memset(di, 0, sizeof(*di));
            di->m = 2;
            di->idx = malloc(sizeof(struct data_block)*di->m);
            di->idx[0].start = start-1;
            di->idx[0].end = end;
            kh_val((kh_idx_t*)idx->dict, k) = di;
            di->n++;
            idx->n++;
            //debug_print("%s\t%d\t%d",s, di->idx[0].start, di->idx[0].end);
            
        }
        else {
            struct data_index *di = kh_val((kh_idx_t*)idx->dict, k);
            if (di->m == di->n) {
                di->m = di->m<<1;
                di->idx = realloc(di->idx, di->m*sizeof(struct data_block));
            }
            di->idx[di->n].start = start-1;
            di->idx[di->n].end = end;
            //debug_print("%s\t%d\t%d",s, di->idx[di->n].start, di->idx[di->n].end);
            di->n++;
        }
        //debug_print("%s\t%d\t%d", s, start, end);
        free(s);
    }
    
    fclose(fp);
    return idx;
}

static void *run_it(void *_d, int i)
{
    struct data_index *idx = (struct data_index*)_d;    
    struct bseq_pool *p = bseq_pool_init();

    FILE *fp = args.fp[i];
    
    int k;   
    for (k = 0; k < idx->n; ++k) {
        struct bseq *b = read_file_block(fp, &idx->idx[k]);
        bseq_pool_push(b, p);
    }
    if (args.dedup == 1)
        bseq_pool_dedup(p);
    return p;        
}

static void write_out(struct bseq_pool *p)
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
    }

    args.out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    return 0;
}

int fsort(int argc, char ** argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();

    struct FILE_tag_index *idx = build_file_index(args.input_fname);

    LOG_print("Build index time: %.3f sec; CPU: %.3f sec", realtime()-t_real, cputime());

    // thread cache
    args.fp = malloc(args.n_thread*sizeof(void*));
    int i;
    for (i = 0; i < args.n_thread; ++i) args.fp[i] = fopen(args.input_fname, "r");
    

    // test_file_index(idx, args.input_fname);

    int n = args.n_thread;
    struct thread_pool *p = thread_pool_init(n);
    struct thread_pool_process *q = thread_pool_process_init(p, n*2, 1);
    struct thread_pool_result *r;
    
    for (i = 0; i < idx->n; ++i) {
        char *key = idx->key[i];
        khint_t k;
        k = kh_get(idx,(kh_idx_t*)idx->dict, key);
        struct data_index *di = kh_val((kh_idx_t*)idx->dict, k);

        int block;

        do {
            block = thread_pool_dispatch2(p, q, run_it, di, 0);
            if ((r = thread_pool_next_result(q))) {
                struct bseq_pool *d = (struct bseq_pool*)r->data;
                write_out(d);
            }
            thread_pool_delete_result(r, 0);
        }
        while (block==-1);
    }
    thread_pool_process_flush(q);

    while ((r = thread_pool_next_result(q))) {
        struct bseq_pool *d = (struct bseq_pool*)r->data;
        write_out(d);
        thread_pool_delete_result(r,0);
    }

    thread_pool_process_destroy(q);
    thread_pool_destroy(p);

    FILE_tag_destroy(idx);
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
