#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "dict.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include <zlib.h>
#include <ctype.h>
#include <sys/stat.h>

static int usage()
{
    fprintf(stderr, "* Sort fastq records by tags and deduplicate.\n");
    fprintf(stderr, "fastq-sort  in.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, " -tag         Tags, such as CB,UR. Order of these tags is sensitive.\n");
    fprintf(stderr, " -dedup       Remove dna copies with same tags. Only keep reads have the best quality.\n");
    fprintf(stderr, " -list        White list for first tag, usually for cell barcodes.\n");
    fprintf(stderr, " -t           Threads.\n");
    fprintf(stderr, " -o           Output fastq.\n");
    fprintf(stderr, " -m           Memory per thread. [1G]\n");
    fprintf(stderr, " -p           Input fastq is smart pairing.\n");
    fprintf(stderr, " -T PREFIX    Write temporary files to PREFIX.nnnn.fq\n");
    fprintf(stderr, " -dropN       Drop if N found in tags.\n");
    fprintf(stderr, "\n");
    return 1;
}

#define MIN_MEM_PER_THREAD  100 // 10M

static struct args {
    const char *input_fname;
    const char *list_fname;
    const char *output_fname;
    const char *prefix;
    int dropN;    
    int n_thread;
    int dedup;
    int smart_pairing;
    int check_list;
    int mem_per_thread;
    struct dict *tags;
    FILE *out;
    FILE *fp_in;
    struct dict *bcodes;
} args = {
    .input_fname = NULL,
    .list_fname = NULL,
    .output_fname = NULL,
    .prefix = "./fsort_temp",
    .dropN = 0,
    .smart_pairing = 0,
    .dedup = 0,
    .mem_per_thread = 1000000000, // 1G
    .tags = NULL,
    .check_list = 0,
    .out = NULL,
    .fp_in = NULL,
    .bcodes = NULL,
};

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *thread = NULL;
    const char *tags = NULL;
    const char *memory = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-list") == 0) var = &args.list_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-prefix") == 0) var = &args.prefix;
        else if (strcmp(a, "-m") == 0) var = &memory;
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

    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (tags == NULL) error("-tag must be set.");

    kstring_t str = {0,0,0};
    kputs(tags, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    
    args.tags = dict_init();    
    for (i = 0; i <n; ++i)
        dict_push(args.tags, str.s+s[i]);

    free(s); free(str.s);
    
    if (thread) args.n_thread = str2int(thread);
    if (args.n_thread < 1) args.n_thread = 1;
    if (memory) args.mem_per_thread = human2int(memory);

    if (args.mem_per_thread < MIN_MEM_PER_THREAD) args.mem_per_thread = MIN_MEM_PER_THREAD;
    
    if (args.list_fname) {
        args.bcodes = dict_init();
        if (dict_read(args.bcodes, args.list_fname))
            error("Barcode list is empty.");
        args.check_list = 1;
    }

    return 0;
}

#define MAX_OPEN_FILE  100

struct fastq_idx {
    int n, m;
    uint32_t *idx;
};

KHASH_MAP_INIT_STR(name, struct fastq_idx)

struct fastq_stream {
    int n, m;
    char **names;
    kh_name_t *dict;
    size_t l_buf;
    uint8_t *buf;
    char *out_fn;
    int reads; // reads in this block
    int dedup; // deduplicated reads in this block
    int pair;
    int fasta;
};
void fastq_stream_destroy(struct fastq_stream *d)
{
    int i;
    for (i = 0; i < d->n; ++i) {
        khint_t k;
        k = kh_get(name, d->dict, d->names[i]);
        free(kh_val(d->dict, k).idx);
        free(d->names[i]);
    }
    free(d->names);
    kh_destroy(name, d->dict);
    free(d->buf);
    free(d);
}
struct fastq_stream_handler {
    FILE *fp;
    int l_buf;
    uint8_t *buf;
    int pair;
    int fasta;
    int mem_per_thread;
};
struct fastq_stream_handler *fastq_stream_handler_init(const char *fn)
{
    struct fastq_stream_handler *f = malloc(sizeof(*f));
    memset(f, 0, sizeof(*f));
    
    f->fp = fopen(fn, "r");
    if (f->fp == NULL) error("%s : %s.", fn, strerror(errno));
    uint8_t c = fgetc(f->fp);
    switch(c) {
        case '@': f->fasta = 0; break;
        case '>': f->fasta = 1; break;
        default: error("Unknown input format. %s", fn); break;
    }

    fseek(f->fp, 0, SEEK_SET);
    
    return f;
}
void fastq_stream_handler_destroy(struct fastq_stream_handler *f)
{
    fclose(f->fp);
    if (f->l_buf) free(f->buf);
}

#define BUF_SEG 1000000 // 1M

struct fastq_stream *fastq_stream_read_block(struct fastq_stream_handler *fastq)
{
    struct fastq_stream *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(*s));

    s->buf = calloc(fastq->mem_per_thread+ fastq->l_buf, 1);
    memcpy(s->buf, fastq->buf, fastq->l_buf);
    s->l_buf = fastq->l_buf;
    if (fastq->l_buf != 0) {
        free(fastq->buf);
        fastq->l_buf = 0;
    }
    size_t ret = fread(s->buf+s->l_buf, 1, fastq->mem_per_thread, fastq->fp);
    s->l_buf += ret;
    
    if (ret == 0 && s->l_buf == 0) {
        free(s);
        return NULL;
    }
    
    if (ret < fastq->mem_per_thread) {
        return s;
    }
    fastq->buf = calloc(BUF_SEG, 1);
    ret = fread(fastq->buf, 1, BUF_SEG, fastq->fp);
    if (ret < BUF_SEG) {
        s->buf = realloc(s->buf, s->l_buf + ret);
        memcpy(s->buf+s->l_buf, fastq->buf, ret);
        s->l_buf+=ret;
        free(fastq->buf);
        fastq->l_buf = 0;
        return s;
    }
    fastq->l_buf = ret;
    int i = 0;
    uint8_t *p = fastq->buf;
    int l; // push to buf
    for ( ;;) {
        for (; i < fastq->l_buf; ++i) if (p[i] == '\n') break;
        ++i; // skip \n        
        if (fastq->fasta == 1) {
            if (p[i] == '>') { // header
                l = i;
                if (l == 0) return s; //
                
                s->buf = realloc(s->buf, s->l_buf + l);
                memcpy(s->buf+s->l_buf, fastq->buf, l);
                s->l_buf += l;
                
                memmove(fastq->buf, fastq->buf+l, fastq->l_buf -l);
                fastq->l_buf -= l;
                return s;
            }
            else continue;
        }
        // fastq format
        if (p[i] == '@') { // this line could be read name or quality string
            l = i;
            for (; i < fastq->l_buf; ++i) if (p[i] == '\n') break;
            i++;
            if (p[i] == '@') { // last line is quality string
                l = i;
                s->buf = realloc(s->buf, s->l_buf + l);
                memcpy(s->buf+s->l_buf, fastq->buf, l);
                s->l_buf += l;
                
                memmove(fastq->buf, fastq->buf+l, fastq->l_buf -l);
                fastq->l_buf -= l;
                return s;                    
            }
            // DNA sequence
            for (; i < fastq->l_buf; ++i) if (p[i] == '\n') break;
            i++;
            if (p[i] != '+') error("Unknown format.");
            
            s->buf = realloc(s->buf, s->l_buf + l);
            memcpy(s->buf+s->l_buf, fastq->buf, l);
            s->l_buf += l;
            
            memmove(fastq->buf, fastq->buf+l, fastq->l_buf -l);
            fastq->l_buf -= l;
            return s;                    
        }
        else continue;
    }
}

static char *query_tags(uint8_t *p, struct dict *dict, int *e)
{
    if (p == NULL) { *e = 1; return NULL; }
    char **val = calloc(dict_size(dict),sizeof(char*));
    int i;
    *e = 0;
    for (i = 0; p[i] != '\n' && p[i] != '\0';) {
        if (p[i] == '\0') { *e = 1; break; }
        if (p[i++] == '|') {
            if (p[i++] == '|') {
                if (p[i++] == '|') {
                    char tag[3];
                    tag[0] = p[i++];
                    tag[1] = p[i++];
                    tag[2] = '\0';
                    if (p[i++] != ':') continue; // not a tag format
                    int idx = dict_query(dict, tag);
                    if (idx == -1) continue;
                    i++; // skip tag character
                    if (p[i++] != ':') val[idx] = strdup("1"); // flag tag
                    else {
                        int j;
                        for (j = i; p[j] != '\n' && p[j] != '\0' && p[j] != '|'; ++j);
                        if (p[j] == '\0') { *e = 1; break; }
                        char *v = malloc(j-i+1);
                        memcpy(v, p+i, j-i);
                        v[j-i] = '\0';
                        val[idx] = v;
                    }
                }
                else if (p[i] == '\0') { *e = 1; break; }
            }
            else if (p[i] == '\0') { *e = 1; break; }
        }
        else if (p[i] == '\0') { *e = 1; break; }
    }

    int empty = 0;
    kstring_t str = {0,0,0};
    for (i = 0; i < dict_size(dict); ++i) {
        if (val[i] == NULL) {
            empty = 1;
            continue;
        }
        if (empty == 0) kputs(val[i], &str);
        free(val[i]);
    }
    free(val);
        
    if (empty == 1) {
        if (str.m) free(str.s);
        return NULL;
    }
    return str.s;
}

int name_cmp(const void *a, const void *b)
{
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    return strcmp(*ia, *ib);
}
char *get_rname(uint8_t *p)
{
    int i;
    for (i = 0; p[i] != '\n' && p[i] != '\0'; ++i);
    char *r = malloc(i+1);
    memcpy(r, p, i);
    r[i] = '\0';
    return r;
}
void sort_block(struct fastq_stream *buf)
{
    buf->dict = kh_init(name);
    uint8_t *p = buf->buf;
    uint8_t *e = buf->buf + buf->l_buf;
    int offset = 0;
    int end;
    for ( ; p != e; ) {
        char *name = query_tags(p, args.tags, &end);
        char *rname = get_rname(p);
        if (name == NULL) error("No tag found at %s", rname);
        
        if (buf->n == buf->m) {
            buf->m = buf->m == 0 ? 1024 : buf->m<<1;
            buf->names = realloc(buf->names, buf->m*sizeof(void*));
        }
        
        khint_t k = kh_get(name, buf->dict, name);
        if (k == kh_end(buf->dict)) {
            buf->names[buf->n] = strdup(name);
            int ret;
            k = kh_put(name, buf->dict, buf->names[buf->n], &ret);
            struct fastq_idx *idx = &kh_val(buf->dict, k);
            memset(idx, 0, sizeof(*idx));
            buf->n++;
        }
        free(name);

        struct fastq_idx *idx = &kh_val(buf->dict, k);
        if (idx->n == idx->m) {
            idx->m = idx->m == 0 ? 10 : idx->m<<1;
            idx->idx = realloc(idx->idx, 4*idx->m);
        }
        idx->idx[idx->n] = offset;
        idx->n++;
        int i = 0;
        for ( ; p + i != e; ++i) if (p[i] == '\n') break; // name
        i++; // skip \n
        for ( ; p + i != e; ++i) if (p[i] == '\n') break; // sequence
        i++;
        if (buf->fasta) {
            p[i-1] = '\0';
        }
        else {
            if (p[i] != '+') {
                error("Format is unrecognised.");
            }
            for ( ; p + i != e; ++i) if (p[i] == '\n') break; // + line
            i++; // skip \n           
            for ( ; p + i != e; ++i) if (p[i] == '\n') break; // qual
            i++;
            if (buf->pair) {
                char *rname2 = get_rname(p);
                if (strcmp(rname, rname2) != 0) error("Inconsistance read name, %s vs %s", rname, rname2);
                free(rname2);
                for ( ; p + i != e; ++i) if (p[i] == '\n') break; // read name
                i++; // skip \n
                for ( ; p + i != e; ++i) if (p[i] == '\n') break; // sequence
                i++;
                for ( ; p + i != e; ++i) if (p[i] == '\n') break; // +
                i++; // skip \n
                for ( ; p + i != e; ++i) if (p[i] == '\n') break; // qual
                i++;
                p[i-1] = '\0';
            }
            else {
                p[i-1] = '\0';
            }
        }
        free(rname);
        offset += i;
        p += i;
        // p = buf->buf + offset;
    }

    qsort(buf->names, buf->n, sizeof(char*), name_cmp);

    // todo: dedup
}

int write_file(struct fastq_stream *buf)
{
    FILE *out = fopen(buf->out_fn, "w");
    if (out == NULL) error("%s : %s.", buf->out_fn, strerror(errno));
    
    int i;
    for (i = 0; i < buf->n; ++i) {
        khint_t k;
        k = kh_get(name, buf->dict, buf->names[i]);
        assert(k != kh_end(buf->dict));
        struct fastq_idx *idx = &kh_val(buf->dict, k);
        int j;
        for (j = 0; j < idx->n; ++j) {
            fputs((char*)(buf->buf+idx->idx[j]), out);
            fputc('\n', out);
        }
    }

    fclose(out);
    // fastq_stream_destroy(buf);
    return 0;
}

struct fastq_stream_reader {
    FILE *fp;
    int blength; // block length
    int n, m; // alloced length of the block
    uint8_t *buf;
    char *name;
    int pair;
    int fasta;
};

struct fastq_stream_reader *fastq_stream_reader_init(const char *fn, int pair)
{
    struct fastq_stream_reader *r = malloc(sizeof(*r));
    memset(r, 0, sizeof(struct fastq_stream_reader));
    r->fp = fopen(fn, "r");
    if (r->fp == NULL) error("%s : %s.", fn, strerror(errno));
    char c = fgetc(r->fp);
    switch(c) {
        case '>': r->fasta = 1; break;
        case '@': r->fasta = 0; break;
        default: error("Unknown file format."); break;
    }

    fseek(r->fp, 0, SEEK_SET);

    r->pair = pair;
    return r;
}

void fastq_stream_reader_close(struct fastq_stream_reader *r)
{
    fclose(r->fp);
    if (r->name) free(r->name);
    if (r->m) free(r->buf);
    free(r);
}
uint8_t *read_next_read(struct fastq_stream_reader *r, uint8_t *p, int *n, int *e)
{
    int i = 0;
    *n = 0;
    for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // rname
    if (p[i] == '\0') goto truncated_buf;
    else ++i; // skip \n

    for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // seq
    if (p[i] == '\0') goto truncated_buf;
    else ++i; // skip \n
    
    if (r->fasta == 1) {
        *n = i;
        return p+i;
    }
    if (p[i] != '+') {
        error("Bad format.");
    }
    
    for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // +\n
    if (p[i] == '\0') goto truncated_buf;
    else ++i; // skip \n

    for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // qual
    if (p[i] == '\0') goto truncated_buf;
    else ++i; // skip \n

    if (r->pair) {
        for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // rname
        if (p[i] == '\0') goto truncated_buf;
        else ++i; // skip \n

        for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // seq
        if (p[i] == '\0') goto truncated_buf;
        else ++i; // skip \n

        if (p[i] != '+') error("Bad format.");

        for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // +\n
        if (p[i] == '\0') goto truncated_buf;
        else ++i; // skip \n

        for ( ; p[i] != '\n' && p[i] != '\0'; ++i); // qual
        if (p[i] == '\0') goto truncated_buf;
        else ++i; // skip \n
    }

    *n = i;
    return p+i;

  truncated_buf:
    *e = 1;
    return NULL;
}

int fastq_stream_reader_sync(struct fastq_stream_reader *r)
{
    if (r->blength != 0) error("Call sync function after print.");
    uint8_t *p = r->buf;
    if (p==NULL) {
      load_buffer:
        r->blength = 0;
        if (r->fp) { // cache more reads
            int new_alloc = r->n + BUF_SEG+1;
            if (new_alloc > r->m) {
                r->buf = realloc(r->buf, new_alloc);
                r->m = new_alloc;
            }
            int ret;
            ret = fread(r->buf+r->n, 1, BUF_SEG, r->fp);
            r->n += ret;
            r->buf[r->n] = '\0';
            if (ret < BUF_SEG) {
                fclose(r->fp);
                r->fp = NULL;
            }
        }
        else if (r->n == 0) return 1;
        p = r->buf;
    }
    
    if (r->name) {
        free(r->name);
        r->name = NULL;
    }
    
    int end = 0;
    r->name = query_tags(p, args.tags, &end);
    if (end) goto load_buffer;

    for ( ;;) {
        if (r->blength == r->n) break;
        int block;
        p = read_next_read(r, p, &block, &end);
        r->blength += block;
        if (end) goto load_buffer;
        if (*p == '\0') {
            r->buf[r->blength-1] = '\0';
            break;
        }
        char *name = query_tags(p, args.tags, &end);        
        if (end) goto load_buffer;

        if (strcmp(r->name, name) != 0) {
            r->buf[r->blength-1] = '\0';
            free(name);
            break;
        }
        free(name);
    }
    return 0;
    //return r->n == 0;
}

int fastq_stream_print(FILE *out, struct fastq_stream_reader *r)
{
    if(r->blength==0) {
        assert(r->n != 0);
        return 1;
    }
    fputs((char*)r->buf, out);
    fputc('\n', out);
    r->n = r->n - r->blength;
    if (r->n == 0) return 1;
    
    memmove(r->buf, r->buf+r->blength, r->n);
    r->blength = 0;
    r->buf[r->n] = '\0';
    free(r->name);
    r->name = NULL;
    return fastq_stream_reader_sync(r);
}
int cmpfunc(const void *a, const void *b)
{
    struct fastq_stream_reader**ia = (struct fastq_stream_reader**)a;
    struct fastq_stream_reader**ib = (struct fastq_stream_reader**)b;
    if (*ia == NULL) return 1;
    if (*ib == NULL) return -1;
    
    return strcmp((*ia)->name,(*ib)->name);
}
int merge_files(const char *fn, int n_file, char **files, int pair)
{
    if (n_file == 1) {
        rename(files[0], fn);
        return 0;
    }
    FILE *out = fn == NULL ? stdout :  fopen(fn, "w");
    if (out == NULL) error("%s : %s.", fn, strerror(errno));
    
    struct fastq_stream_reader **readers = malloc(n_file*sizeof(void *));
    int i;
    for (i = 0; i < n_file; ++i) {
        readers[i] = fastq_stream_reader_init(files[i], pair);
        fastq_stream_reader_sync(readers[i]);
    }

    for (;;) {
        qsort(readers, n_file, sizeof(void*), cmpfunc);
        int j;
        int finished = 1;
        char *last_name = NULL;
        for (j = 0; j < n_file; ++j) {
            if (readers[j] == NULL) continue;
            finished = 0;
            assert(readers[j]->name != NULL);

            if (last_name == NULL) {    
                last_name = strdup(readers[j]->name);
                int ret;
                ret = fastq_stream_print(out, readers[j]);
                if (ret == 1) {
                    fastq_stream_reader_close(readers[j]);
                    readers[j] = NULL;
                }
            }
            else if (strcmp(last_name, readers[j]->name) == 0) {
                int ret;
                ret = fastq_stream_print(out, readers[j]);
                if (ret == 1) {
                    fastq_stream_reader_close(readers[j]);
                    readers[j] = NULL;
                }
            }
            else {
                free(last_name);
                break;
            }
        }
        if (finished == 1) break;
    }
    fclose(out);
    return 0;
}

static void *run_it(void *d)
{
    sort_block((struct fastq_stream *)d);
    return d;
}

int fsort(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();
    // init fastq handler
    struct fastq_stream_handler *fastq = fastq_stream_handler_init(args.input_fname);
    fastq->pair = args.smart_pairing;
    fastq->mem_per_thread = args.mem_per_thread;
    
    int n_file = 0;
    int m_file = 0;
    char **files = NULL;
    
    hts_tpool *pool = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(pool, args.n_thread*2, 0);
    hts_tpool_result *r;

    for ( ;;) {

        struct fastq_stream *stream = fastq_stream_read_block(fastq);
        if (stream == NULL) break;

        if (n_file == m_file) {
            m_file = m_file == 0 ? 10 : m_file<<1;
            files = realloc(files, m_file*sizeof(char*));
        }
        files[n_file] = calloc(strlen(args.prefix)+20,1);
        sprintf(files[n_file], "%s.%.4d.fq", args.prefix, n_file);
        stream->out_fn = files[n_file];
        n_file++;
        
        int block;
        do {
            block = hts_tpool_dispatch2(pool, q, run_it, stream, 1);
            if ((r = hts_tpool_next_result(q))) {    
                struct fastq_stream *d = (struct fastq_stream*)hts_tpool_result_data(r);
                write_file(d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);
        
    while ((r = hts_tpool_next_result(q))) {
        struct fastq_stream *d = (struct fastq_stream*)hts_tpool_result_data(r);
        write_file(d);
        
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(pool);

    
    merge_files(args.output_fname, n_file, files, args.smart_pairing);

    int i;
    for (i = 0; i < n_file; ++i) {
        unlink(files[i]);
        free(files[i]);
    }
    free(files);
    fastq_stream_handler_destroy(fastq);

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());

    return 0;
}
