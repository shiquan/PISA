#include "utils.h"
#include "number.h"
#include "fml.h"
#include "fastq.h"
#include "htslib/thread_pool.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

int usage()
{
    fprintf(stderr, "assem  in.fq\n");
    fprintf(stderr, "  -t       Threads.\n");
    fprintf(stderr, "  -tag     Tags to specify each block.\n");
    fprintf(stderr, "  -p       Smart pairing.\n");
    fprintf(stderr, "  -i       Index file of in.fq\n");
    return 1;
}
struct args {
    const char *input_fname;
    const char *output_fname;
    int n_tag;
    char **tags;
    int smart_pairing;
    
    kseq_t *ks;
    char *last_name;
    gzFile fp_in;
    FILE *out;

    fml_opt_t *assem_opt;
    int n_thread;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .n_tag = 0,
    .tags = NULL,
    .smart_pairing = 0,
    .ks = NULL,
    .last_name = NULL,
    .fp_in = NULL,
    .out = NULL,
    .assem_opt = NULL,
    .n_thread = 1,
};
void memory_release()
{
    int i;
    for (i = 0; i < args.n_tag; ++i) free(args.tags[i]);
    free(args.tags);
    kseq_destroy(args.ks);
    gzclose(args.fp_in);
    fclose(args.out);
}
static void assem_opt_init(fml_opt_t *opt)
{
    fml_opt_init(opt);
    opt->min_asm_ovlp = 10;
    opt->mag_opt.flag = MAG_F_NO_SIMPL | MAG_F_AGGRESSIVE;
};

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
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
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

    args.out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    args.assem_opt = malloc(sizeof(fml_opt_t));
    assem_opt_init(args.assem_opt);

    args.fp_in = gzopen(args.input_fname, "r");
    CHECK_EMPTY(args.fp_in, "%s : %s.", args.input_fname, strerror(errno));

    args.ks = kseq_init(args.fp_in);
                
    return 0;
}

struct read_block {
    char *name;
    bseq1_t *b;
    int n, m;
};
void read_block_destory(struct read_block *r)
{
    free(r->name);
    int i;
    for (i = 0; i < r->n; ++i) {
        free(r->b[i].seq);
        free(r->b[i].qual);
    }
    free(r->b);
    free(r);    
}
static char *generate_names(char **names)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < args.n_tag; ++i) kputs(names[i], &str);
    for (i = 0; i < args.n_tag; ++i) {
        kputs("|||", &str);
        kputs(args.tags[i], &str);
        kputs("|||", &str);
        kputs(names[i], &str);
    }
    return str.s;
}
struct read_block *read_block()
{
    struct read_block *b =  malloc(sizeof(*b));
    memset(b, 0, sizeof(*b));
    b->m = 2;
    b->b = malloc(sizeof(bseq1_t)*b->m);

    if (args.last_name) {
        b->name = strdup(args.last_name);

        bseq1_t *bb = &b->b[0];
        bb->seq = strdup(args.ks->seq.s);
        bb->qual =  strdup(args.ks->qual.s);
        bb->l_seq = args.ks->seq.l;
        b->n++;
        free(args.last_name);
        args.last_name = NULL;
    }

    // read name:  TAG1TAGS|||TAG1|||VAL1|||TAG2|||VAL2
    while(kseq_read(args.ks)>=0) {
        if (b->n == b->m) {
            b->m = b->m<<1;
            b->b = realloc(b->b, sizeof(bseq1_t)*b->m);
        }

        char **name = fastq_name_pick_tags(args.ks->name.s, args.n_tag, args.tags);
        char *n = generate_names(name);
        int i;
        for (i = 0; i <args.n_tag; ++i) free(name[i]);
        free(name);
        if (b->name == NULL) b->name = strdup(n);
        // if (args.last_name == NULL) args.last_name = strdup(n); // first line
        if (strcmp(n, b->name) == 0) {
            bseq1_t *b1 = &b->b[b->n++];
            b1->seq = strdup(args.ks->seq.s);
            b1->qual = strdup(args.ks->qual.s);
            b1->l_seq = args.ks->seq.l;
            free(n);
        }
        else {
            //if (args.last_name) free(args.last_name);
            args.last_name = n;
            break;
        }
    }
    if (b->n == 0) {
        free(b); return NULL;
    }
    return b;
}

static char *rend_bseq(struct read_block *b)
{
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < b->n; ++i) {
        bseq1_t *b1 = &b->b[i];
        kputc('@', &str);
        kputw(i, &str);kputc('_', &str);
        kputs(b->name, &str);
        kputc('\n', &str);
        kputs(b1->seq, &str); kputc('\n', &str);
        kputc('+', &str);  kputc('\n', &str);
        kputs(b1->qual, &str); kputc('\n', &str);
    }
    read_block_destory(b);
    kputs("", &str);
    return str.s;
}

static char *rend_utg(char *name, fml_utg_t *utg, int n)
{
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < n; ++i) {
        fml_utg_t *u = &utg[i];
        kputc('@', &str);
        kputw(i, &str);kputc('_', &str);
        kputs(name, &str);
        kputc('\t', &str);
        kputs("assembled",&str);
        kputc('\n', &str);
        kputsn(u->seq, u->len, &str); kputc('\n', &str);
        kputc('+', &str);  kputc('\n', &str);
        kputsn(u->cov, u->len, &str); kputc('\n', &str);
    }
    kputs("", &str);
    return str.s;
}

static void *run_it(void *_d)
{
    char *out = NULL;
    struct read_block *b = (struct read_block*)_d;
    //if (b->n == 1) {
    //out = rend_bseq(b);
//}
    //else {
    if (b->n > 1) {
        int n_utg;
        fml_utg_t *utg = fml_assemble(args.assem_opt, b->n, b->b, &n_utg);
        out = rend_utg(b->name, utg, n_utg);
        fml_utg_destroy(n_utg, utg);
    }
    //read_block_destory(b);
    return out;
}
static void write_out(char *s)
{
    fputs(s, args.out);
    free(s);
}
int assem(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        //int nseq;
        struct read_block *b = read_block();
        if (b == NULL) break;
        if (b->n == 1) continue;
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                char *s = (char*)hts_tpool_result_data(r);
                write_out(s);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }

    hts_tpool_process_flush(q);
    
    while ((r = hts_tpool_next_result(q))) {
        char *s = (char*)hts_tpool_result_data(r);
        write_out(s);
        hts_tpool_delete_result(r, 0);
    }
        
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    return 0;
}

int main(int argc, char **argv)
{
    return assem(argc, argv);
}
