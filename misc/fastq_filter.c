#include "utils.h"
#include "fastq.h"
#include "number.h"
#include <string.h>
#include "htslib/thread_pool.h"

static int usage()
{
    fprintf(stderr, "fastq_filter in1.fq [in2.fq]\n");
    fprintf(stderr, " -p         Input is smart pairing.\n");
    fprintf(stderr, " -k         Apply BGISEQ filter rule.\n");
    fprintf(stderr, " -q [30]    Filter if average quality smaller than this value.\n");
    fprintf(stderr, " -t [5]     Threads.\n");
    return 1;
}
static struct args {
    const char *r1_fname;
    const char *r2_fname;
    
    int smart_pairing;
    int bgi_filter;
    int qual;
    int n_thread;

    struct fastq_handler *fq;

    uint64_t raw_reads;
    uint64_t pass_reads;
} args = {
    .r1_fname = NULL,
    .r2_fname = NULL,
    .smart_pairing = 0,
    .bgi_filter = 0,
    .qual = 30,
    .n_thread = 5,
    .fq = NULL,
    .raw_reads = 0,
    .pass_reads = 0,
};
    
static int parse_args(int argc, char **argv)
{
    if (argc == 1) return usage();
    int i;
    const char *qual = NULL;
    const char *thread = NULL;
    for (i = 1; i <argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-p") == 0) {
            args.smart_pairing = 1;
            continue;
        }
        else if (strcmp(a, "-k") == 0) {
            args.bgi_filter = 1;
            continue;
        }
        else if (strcmp(a, "-q") == 0) var = &qual;
        else if (strcmp(a, "-t") == 0) var = &thread;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.r1_fname == NULL) {
            args.r1_fname = a;
            continue;
        }
        if (args.r2_fname == NULL) {
            args.r2_fname = a;
            continue;
        }
        error("Unknown argument, %s.",a);
    }

    CHECK_EMPTY(args.r1_fname, "No input fastq");

    if (qual) args.qual = str2int((char*)qual);
    if (args.qual < 0) args.qual = 0;

    if (thread) args.n_thread = str2int((char*)thread);
    if (args.n_thread < 1) args.n_thread = 1;

    args.fq = fastq_handler_init(args.r1_fname, args.r2_fname, args.smart_pairing, 80000000);
    CHECK_EMPTY(args.fq, "Failed to load fastq.");
    
    return 0;
}

#define FQ_PASS 0
#define FQ_FAIL 1

void *run_it(void *_d)
{
    struct bseq_pool *p = (struct bseq_pool*)_d;
    struct args *opt = (struct args*)p->opts;

    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        b->flag = FQ_PASS;
        if (opt->bgi_filter) {
            int k;
            int bad_bases = 0;
            for (k = 0; k < 15 && k < b->l0; ++k) 
                if (b->q0[k]-33<10) bad_bases++;
                
            if (bad_bases > 2) {
                b->flag = FQ_FAIL;
                continue;
            }
            if (b->l1 >0) {
                for (k = 0; k < 15 && k < b->l1; ++k) 
                    if (b->q1[k]-33<10) bad_bases++;
                if (bad_bases >2) {
                    b->flag = FQ_FAIL;
                    continue;
                }
            }
        }
        if (opt->qual > 0 && b->flag == FQ_PASS) {
            int k;
            int ave = 0;
            for (k = 0; k < b->l0; k++) ave += b->q0[k]-33;
            if (ave/k < opt->qual) {
                b->flag = FQ_FAIL;
                continue;
            }
            if (b->l1 > 0) {
                ave = 0;
                for (k = 0; k < b->l1; k++) ave += b->q1[k]-33;
                if (ave/k < opt->qual) {
                    b->flag = FQ_FAIL;
                    continue;
                }
            }
        }
    }
    return p;
}
static void write_out(void *_d)
{
    struct bseq_pool *p = (struct bseq_pool*)_d;
    struct args *opt = (struct args*)p->opts;
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        opt->raw_reads++;
        if (b->flag == FQ_PASS) {
            opt->pass_reads++;
            printf("@%s\n%s\n+\n%s\n", b->n0, b->s0, b->q0);
            if (b->l1 > 0) printf("@%s\n%s\n+\n%s\n", b->n0, b->s1, b->q1);
        }
    }
    bseq_pool_destroy(p);
}
static void memory_release()
{
    fastq_handler_destory(args.fq);
}
int main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    for (;;) {
        struct bseq_pool *b = bseq_read(args.fq, &args);
        if (b == NULL) break;
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 0);
            if ((r = hts_tpool_next_result(q))) {                
                struct bseq_pool *d = (struct bseq_pool*)hts_tpool_result_data(r);
                write_out(d);
            }
            hts_tpool_delete_result(r, 1);            
        }
        while (block == -1);
    }
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct bseq_pool *d = (struct bseq_pool*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    memory_release();

    fprintf(stderr, "Raw fragment : %"PRIu64"\n", args.raw_reads);
    fprintf(stderr, "QC passed : %"PRIu64"\n", args.pass_reads);
    return 0;
}
