#include "utils.h"
#include "number.h"
#include "thread_pool.h"
#include "fastq.h"
#include <string.h>
#include "htslib/kstring.h"

static int usage()
{
    fprintf(stderr, "* This program used to trim the mosic ends and merge paired reads for scLFR library.\n");
    fprintf(stderr, "Usage: LFR_cleanup input.fq\n");
    return 1;
}

static struct args {
    const char *input_fname;
    const char *output_fname;
    int n_thread;
    FILE *out;
    uint8_t *me_enc;
    uint8_t *rev_enc;
    struct fastq_handler *fastq;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .n_thread = 1,
    .out = NULL,
    .me_enc = NULL,
    .rev_enc = NULL,
    .fastq = NULL,
};

static void memory_release()
{
    fclose(args.out);
    free(args.me_enc);
    fastq_handler_destory(args.fastq);
}
extern int check_compl(const char *s1, const char *s2, int l);
extern uint8_t *nt4_enc(char *s, int l);
extern int check_overlap(const int lr, const int lq, char const *r, uint8_t const *qry);
extern int merge_paired(const int l_seq1, const int l_seq2, char const *s1, char const *s2, char const *q1, char const *q2, char **_seq, char **_qual);
extern uint8_t *rev_enc(uint8_t *s, int l);
static int parse_args(int argc, char **argv)
{
    int i;
    const char *thread = NULL;
    for (i = 1; i < argc;) {

        const char *a = argv[i++];
        const char **var = 0;
        
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return usage();
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1]) error("Unknown parameter. %s", a);
        
        if (args.input_fname == 0) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument: %s, use -h see help information.", a);
    }

    if (thread) args.n_thread = str2int((char*)thread);
    if (args.input_fname == NULL && (!isatty(fileno(stdin))))
        args.input_fname = "-";    
    if (args.input_fname == NULL) error("Fastq file(s) must be set!");

    args.fastq = fastq_handler_init(args.input_fname, NULL, 1, 10000);

    if (args.output_fname != NULL) {
        args.out = fopen(args.output_fname, "w");
        if (args.out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    }
    else args.out = stdout;
    const char *ME = "CTGTCTCTTATACACATCT";
    args.me_enc = nt4_enc(ME,19);
    args.rev_enc = rev_enc(args.me_enc, 19);
    return 0;
}

static void write_out(void *_data)
{
    char *s = (char *)_data;
    if (s) {
        fputs(s, args.out);
        free(s);
    }
}
char *trim_ends(struct bseq_pool *p)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        int l;
        // all reverse ME sequence should not be detected
        l = check_overlap(b->l0, 19, b->s0, args.rev_enc);
        if (l != -1) continue;
        l = check_overlap(b->l1, 19, b->s1, args.rev_enc);
        if (l != -1) continue;
        
        l = check_overlap(b->l0, 19, b->s0, args.me_enc);        
        if (l == -1) {
            char *s, *q;
            if (merge_paired(b->l0, b->l1, b->s0, b->s1, b->q0, b->q1, &s, &q) == 0) {
                kputc('@', &str);
                kputs(b->n0, &str); kputc('\n', &str);
                kputs(s, &str); kputs("\n+\n", &str);
                kputs(q, &str); kputc('\n', &str);
                free(s); free(q);
            }
            else {
                kputc('@', &str);
                kputs(b->n0, &str); kputc('\n', &str);
                kputs(b->s0, &str); kputs("\n+\n", &str);
                kputs(b->q0, &str); kputc('\n', &str);
                kputc('@', &str);
                kputs(b->n0, &str); kputc('\n', &str);
                kputs(b->s1, &str); kputs("\n+\n", &str);
                kputs(b->q1, &str); kputc('\n', &str);
            }
        }
        else {
            if (l < 40) continue; // to short;
            //char *s0 = b->s0;
            // char *s1 = b->s1 + (b->l1 - l);
            if (check_compl(b->s0, b->s1, l) == 0) {                                
                kputc('@', &str);
                kputs(b->n0, &str); kputc('\n', &str);
                kputsn(b->s0, l, &str); kputs("\n+\n", &str);
                kputsn(b->q0, l, &str); kputc('\n', &str);
            }
        }
    }

    bseq_pool_destroy(p);
    return str.s;
}

static void *run_it(void *_p, int idx)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    return trim_ends(p);
}

int LFR_cleanup(int argc, char **argv)
{
    if (parse_args(argc, argv)) return 1;
    
    double t_real;
    t_real = realtime();

    void *opts = &args;
    int nt = args.n_thread;
    
    struct thread_pool *p = thread_pool_init(nt);
    struct thread_pool_process *q = thread_pool_process_init(p, nt*2, 0);
    struct thread_pool_result *r;

    for (;;) {
        struct bseq_pool *b = fastq_read(args.fastq, &args);
        if (b == NULL) break;
        int block;
        do {
            block = thread_pool_dispatch2(p, q, run_it, b, 1);
            if ((r = thread_pool_next_result(q))) {
                struct bseq_pool *d = (struct bseq_pool *)r->data;
                write_out(d);
            }
            thread_pool_delete_result(r, 0);
        }
        while (block == -1);
    }
    thread_pool_process_flush(q);

    while ((r = thread_pool_next_result(q))) {
        struct bseq_pool *d = (struct bseq_pool*)r->data;
        write_out(d);
        thread_pool_delete_result(r, 0);
    }

    thread_pool_process_destroy(q);
    thread_pool_destroy(p);
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    
    return 0;
}


int main(int argc, char **argv)
{
    return LFR_cleanup(argc, argv);
}

