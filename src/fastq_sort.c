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
    // fprintf(stderr, " -dropN       Drop if N found in tags.\n");
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
    int paired;
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
    .prefix = NULL,
    .dropN = 0,
    .paired = 0,
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
            args.paired = 1;
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
    if (args.output_fname == NULL) error("Output file should be specified.");
    
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


    if (args.prefix == NULL)
        args.prefix= args.output_fname;
        
    return 0;
}

struct record_offset {
    int n, m;
    uint64_t *offsets;
};

struct read_block {
    struct dict *dict;
    int n, m;
    struct record_offset *idx;
    int max;
    uint8_t *data;
};

struct fastq_idx {
    int n;
    char **name;
    long int *offset;
};

void fastq_idx_destroy(struct fastq_idx *i)
{
    int k;
    for (k = 0; k < i->n; ++k) free(i->name[k]);
    free(i->name);
    free(i->offset);
    free(i);
}
    
struct fastq_node {
    char *fn;
    FILE *fp;
    struct fastq_idx *idx;
    int i;
    char *name; // point to idx::name[i]
    int n, m;
    char *buf;
};

void fastq_node_clean(struct fastq_node *n)
{
    free(n->fn);
    fclose(n->fp);
    fastq_idx_destroy(n->idx);
    if (n->m) free(n->buf);
    memset(n, 0, sizeof(*n));
    //free(n);
    // n = NULL;
}

struct fastq_stream {
    struct read_block *r;
    struct fastq_node *n;
};

struct read_block *read_block_file(FILE *fp, int max, int paired)
{
    struct read_block *r = malloc(sizeof(*r));
    memset(r, 0, sizeof(*r));
    int fasta;

    kstring_t str = {0,0,0};
    
    int record = 0;
    for (;;) {
        char c = fgetc(fp);
        if (c == EOF) break;
        if (c == '@') fasta= 0;
        else if (c == '>') fasta = 1;
        else error("Unknown input format.");

        kputc(c, &str);
        for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); // read name;
        if (c == EOF) error("Truncated file?");
    
        kputc(c, &str); // push \n to buf

        for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); //
        kputc(c, &str); // push \n to buf
    
        if (fasta == 0) {
            c = fgetc(fp);
            if (c != '+') error("Unknown format? %c",c);
            kputc(c, &str); // push + to buf
            for (c = fgetc(fp); c && c!= '\n'; c = fgetc(fp)) kputc(c, &str); // qual name
            kputc(c, &str); // push \n to buf
            for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); //
            kputc(c, &str);
            
            if (paired) {
                // read four more lines
                for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); //
                kputc(c, &str); 
                for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); //
                kputc(c, &str);
                for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); //
                kputc(c, &str);
                for (c = fgetc(fp); c && c != '\n'; c = fgetc(fp)) kputc(c, &str); //
                kputc(c, &str);
            }        
        }
        record++;
        if (str.s[str.l-1] == '\n') str.s[str.l-1] = '\0';
        if (str.l >= max) break;
    }
    if (str.l == 0) {
        free(r);
        return NULL;
    }
    r->data = (uint8_t*) str.s;
    r->max = str.l;

    LOG_print("Read %d records", record);
    return r;
}

static int name_cmp(const void *a, const void *b)
{
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    return strcmp(*ia, *ib);
}
static char *query_tags(uint8_t *p, struct dict *dict)
{
    uint8_t *pp = p;
    for ( ; *pp != '\0' && *pp != '\n'; pp++) {}
    if (*pp == '\0') error("Truncated name.");
    char **val = calloc(dict_size(dict),sizeof(char*));
    int i;
    for (i = 0; p[i] != '\n';) {
        if (p[i++] == '|' && p[i++] == '|' && p[i++] == '|') {
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
                // if (p[j] == '\0') { *e = 1; break; }
                char *v = malloc(j-i+1);
                memcpy(v, p+i, j-i);
                v[j-i] = '\0';
                val[idx] = v;
            }                    
        }
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

void read_block_sort_by_name(struct read_block *r, struct dict *tag)
{
    r->dict = dict_init();
    
    int i;
    for (i = 0; i < r->max; ) {
        char *name = query_tags(r->data+i, tag);
        if (name == NULL) error("Failed to parse names.");
        int id = dict_push(r->dict, name);
        if (id >= r->m) {
            r->m = id+1;
            r->idx = realloc(r->idx, sizeof(struct record_offset)*r->m);
            for (; r->n < r->m; ++r->n) memset(&r->idx[r->n], 0, sizeof(struct record_offset));
            r->n = r->m;
        }

        struct record_offset *off = &r->idx[id];
        if (off->n == off->m) {
            off->m = off->m == 0 ? 4 : off->m*2;
            off->offsets = realloc(off->offsets,off->m*sizeof(uint64_t));
        }
        off->offsets[off->n++] = i;
        for (;r->data[i] != '\0'; ++i) {}
        ++i; // skip \0
    }
    
    char **names = dict_names(r->dict);
    qsort(names, dict_size(r->dict), sizeof(char*), name_cmp);   
}

void read_block_destroy(struct read_block *r)
{
    int i;
    for (i = 0; i < r->n; ++i)
        if (r->idx[i].m) free(r->idx[i].offsets);
    free(r->idx);
    dict_destroy(r->dict);
    free(r->data);
    free(r);
}

struct fastq_idx *write_block_with_idx(const char *fn, struct read_block *r)
{
    FILE *fp = fopen(fn, "w");
    if (fp == NULL) error("%s : %s.", fn, strerror(errno));
    struct fastq_idx *idx = malloc(sizeof(*idx));
    memset(idx, 0, sizeof(*idx));
    idx->n = dict_size(r->dict);
    idx->name = malloc(idx->n*sizeof(char*));
    idx->offset = malloc(idx->n*sizeof(long int));
    
    int i;
    for (i = 0; i < idx->n; ++i) {
        char *name = dict_name(r->dict, i);
        int old_idx = dict_query(r->dict, name);
        struct record_offset *off = &r->idx[old_idx];
        int j;
        for (j = 0; j < off->n; ++j) {
            fputs((char*)(r->data+off->offsets[j]),fp);
            fputc('\n', fp);
        }
        idx->name[i] = strdup(name);
        idx->offset[i] = ftell(fp);
    }
    fclose(fp);
    return idx;
}
static int merge_cmp(const void *a, const void *b)
{
    const struct fastq_node *ia = *(const struct fastq_node **)a;
    const struct fastq_node *ib = *(const struct fastq_node **)b;

    if (ia->name == NULL) return 1;
    if (ib->name == NULL) return -1;

    return strcmp(ia->name, ib->name);
}
int fastq_merge_core(struct fastq_node **node, int n_node, FILE *fp)
{
    int n = n_node;
    char *name = node[0]->name;
    int i;
    for (i = 0; i < n; ++i) {
        struct fastq_node *d = node[i];

        if (d == NULL || d->name == NULL) {
            n=i; // emit all empty handler
            break;
        }
            
        if (strcmp(name, d->name) == 0) {
            d->buf[d->n] = '\0';
            fputs(d->buf, fp);
            // cache next record
            d->i ++;
            if (d->i >= d->idx->n) { // close handler
                unlink(d->fn);
                LOG_print("Unlink %s", d->fn);
                fastq_node_clean(d);
                // debug_print("%p", d);
            }
            else {
                int l;
                l = d->idx->offset[d->i] - d->idx->offset[d->i-1];
                if (l >= d->m) {
                    d->m = l+1;
                    d->buf = realloc(d->buf, d->m);
                }
                d->n = l;
                int ret;
                ret = fread(d->buf, 1, d->n, d->fp);
                d->buf[d->n] = '\0';
                assert(ret == d->n);
                d->name = d->idx->name[d->i];
            }
        }
        else break;
    }
    return n;
}
struct fastq_idx *fastq_merge(struct fastq_node **node, int n_node, const char *fn)
{
    // init
    FILE *fp = fopen(fn, "w");
    if (fp == NULL) error("%s : %s.", fn, strerror(errno));
    int i;
    for (i = 0; i < n_node; ++i) {
        struct fastq_node *d = node[i];
        d->fp = fopen(d->fn, "r");
        // debug_print("%s", d->fn);
        if (d->fp == NULL) error("%s : %s.", d->fn, strerror(errno));
        d->name = d->idx->name[0];
        d->m = d->idx->offset[0] + 1;
        d->buf = malloc(d->m);
        d->n = d->idx->offset[0];
        fread(d->buf, 1, d->n, d->fp);
        d->buf[d->n] = '\0';
    }
    // merge
    int n = n_node;
    int m_idx = 0;
    struct fastq_idx *idx = malloc(sizeof(*idx));
    memset(idx, 0, sizeof(*idx));
    
    for (;;) {
        if (n > 1) 
            qsort(node, n, sizeof(struct fastq_node*), merge_cmp);
        
        if (node[0]->name == NULL) break;

        if (idx->n == m_idx) {
            m_idx = m_idx == 0 ? 1024 : m_idx*2;
            idx->name = realloc(idx->name, m_idx*sizeof(char*));
            idx->offset = realloc(idx->offset, m_idx*sizeof(long int));
        }
        idx->name[idx->n] = strdup(node[0]->name);
        
        n = fastq_merge_core(node, n, fp);
        idx->offset[idx->n] = ftell(fp);
        idx->n++;
    }
    for (i = 0; i < n_node; ++i) free(node[i]);
    fclose(fp);
    LOG_print("Create %s from %d files.", fn, n_node);
    return idx;
}

struct fastq_idx *merge_files(struct fastq_stream *fastqs, int n, const char *fn)
{
    struct fastq_node **nodes = malloc(n*sizeof(struct fastq_node*));
    int i;
    for (i = 0; i < n; ++i) nodes[i] = fastqs[i].n;
    struct fastq_idx *idx = fastq_merge(nodes, n, fn);
    free(nodes);
    return idx;
}

static void *run_it(void *d)
{
    struct fastq_stream *fastq = (struct fastq_stream*)d;
    struct read_block *r = fastq->r;
    struct fastq_node *n = fastq->n;
    read_block_sort_by_name(r, args.tags);
    n->idx = write_block_with_idx(n->fn, r);

    return fastq;
}

static const int max_file_open = 100;

int fsort(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();
    
    FILE *fp = fopen(args.input_fname, "r");
    if (fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    
    int n_file = 0;
    int i_name = 0;
    struct fastq_stream *fastqs = malloc(max_file_open*sizeof(struct fastq_stream));
    
    hts_tpool *pool = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(pool, args.n_thread*2, 0);
    hts_tpool_result *r;

    for ( ;;) {


        if (n_file >= 100) {
            char *name = calloc(strlen(args.prefix)+20,1);
            sprintf(name, "%s.%.4d.fq", args.prefix, i_name);
            i_name++;

            // finished all files
            hts_tpool_process_flush(q);
                   
            while ((r = hts_tpool_next_result(q))) {                
                hts_tpool_delete_result(r, 0);
            }

            struct fastq_idx *idx = merge_files(fastqs, n_file, name);
            memset(fastqs, 0, sizeof(struct fastq_stream)*max_file_open);
            fastqs[0].n = malloc(sizeof(struct fastq_node));
            fastqs[0].n->idx = idx;
            fastqs[0].n->fn = name;
            n_file = 1; // merged file will put as the first record in the list
        }

        struct read_block *b = read_block_file(fp, args.mem_per_thread, args.paired);
        if (b == NULL) break;
        struct fastq_stream *stream = &fastqs[n_file];
        n_file++;

        char *name = calloc(strlen(args.prefix)+20,1);
        sprintf(name, "%s.%.4d.fq", args.prefix, i_name);
        i_name++;

        stream->r = b;
        stream->n = malloc(sizeof(struct fastq_node));
        memset(stream->n, 0, sizeof(struct fastq_node));
        stream->n->fn = name;
        
        int block;
        do {
            block = hts_tpool_dispatch2(pool, q, run_it, stream, 1);
            if ((r = hts_tpool_next_result(q))) {    
                // struct fastq_stream *d = (struct fastq_stream*)hts_tpool_result_data(r);
                // write_file(d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);
        
    while ((r = hts_tpool_next_result(q))) {
        // struct fastq_stream *d = (struct fastq_stream*)hts_tpool_result_data(r);
        // write_file(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(pool);

    char *name = strdup(args.output_fname);
    struct fastq_idx *idx = merge_files(fastqs, n_file, name);
    fastq_idx_destroy(idx);    
    
    // dedup
    
    free(fastqs);


    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());

    return 0;
}
