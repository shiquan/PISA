// convert BGISeq fastq to 10x Input fastq
#include "utils.h"
#include "fastq.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "number.h"

int usage()
{
    // default output read-I1_si-TTTCATGA_lane-008-chunk-001.fastq.gz
    fprintf(stderr, "BGISeqTo10X  -1 read_1.fq.gz -2 read_2.fq.gz -o outdir\n");
    return 1;
}

struct args {
    const char *r1_fname;
    const char *r2_fname;
    const char *output_dir;

    int n_thread;
    
} args = {
    .r1_fname = NULL,
    .r2_fname = NULL,
    .output_dir = NULL,
    .n_thread = 5,
};
int parse_args(int argc, char **argv)
{
    if (argc == 1 ) return usage();
    int i;
    const char **a = NULL;
    const char *thread = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return usage();
        if (strcmp(a, "-1") == 0) var = &args.r1_fname;
        else if (strcmp(a, "-2") == 0) var = &args.r2_fname;
        else if (strcmp(a, "-o") == 0) var = &args.output_dir;
        else if (strcmp(a, "-t") == 0) var = &thread;
                        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        error("Unknown argument: %s, use -h see help information.", a);
    }
    CHECK_EMPTY(args.r1_fname, "Failed to set read 1.");
    CHECK_EMPTY(args.r2_fname, "Failed to set read 2.");
    if (thread) args.n_thread = str2int((char*)thread);
    
    return 0;
}
int main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return 1;
    
    struct fastq_handler *fastq = fastq_handler_init(args.r1_fname, args.r2_fname, 0, 10000);
    kstring_t str1 = {0,0,0};
    kstring_t str2 = {0,0,0};
    if (args.output_dir) {
        kputs(args.output_dir, &str1);
        if (str1.s[str1.l-1] != '/') kputc('/', &str1);
        
        kputs(args.output_dir, &str2);
        if (str2.s[str2.l-1] != '/') kputc('/', &str2);
    }
    kputs("read-RA_si-TTTCATGA_lane-008-chunk-001.fastq.gz", &str1);
    kputs("read-I1_si-TTTCATGA_lane-008-chunk-001.fastq.gz", &str2);
    
    BGZF *fp1 = bgzf_open(str1.s, "w");
    BGZF *fp2 = bgzf_open(str2.s, "w");

    CHECK_EMPTY(fp1, "%s : %s.", str1.s, strerror(errno));
    CHECK_EMPTY(fp2, "%s : %s.", str2.s, strerror(errno));
    
    free(str1.s);
    free(str2.s);
    
    // bgzf_mt(fp1, args.n_thread, 256);
    // bgzf_mt(fp2, args.n_thread, 256);
    
    do {
        struct bseq_pool *p = bseq_read(fastq, &args);
        if ( p == NULL ) break;
        int i;
        kstring_t r1 = {0,0,0};
        kstring_t r2 = {0,0,0};
        
        for (i = 0; i < p->n; ++i) {
            struct bseq *b = &p->s[i];
            kputc('@', &r1);
            kputs(b->n0, &r1);
            kputs(" 1:N:0:0\n", &r1);
            kputs(b->s0, &r1);
            kputs("\n+\n", &r1);
            kputs(b->q0, &r1); 
            kputc('\n', &r1);

            kputc('@', &r1);
            kputs(b->n0, &r1);
            kputs(" 3:N:0:0\n", &r1);
            kputsn(b->s1, b->l1 -8, &r1);
            kputs("\n+\n", &r1);
            kputsn(b->q1, b->l1 -8, &r1);
            kputc('\n', &r1);

            kputc('@', &r2);
            kputs(b->n0, &r2);
            kputs(" 2:N:0:0\n", &r2);
            kputs("TTTCATGA\n+\nAAAAFJ7F\n", &r2); // fack sample name
            bgzf_write(fp1, r1.s, r1.l);
            bgzf_write(fp2, r2.s, r2.l);
            r1.l = 0;
            r2.l = 0;
        }
        bseq_pool_destroy(p);
        free(r1.s);
        free(r2.s);
    } while(1);
    bgzf_close(fp1);
    bgzf_close(fp2);
    fastq_handler_destory(fastq);
    return 0;
}

