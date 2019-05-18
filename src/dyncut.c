#include "utils.h"
#include "number.h"
#include "thread_pool.h"
#include "fastq.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"

#define TRIMME    0
#define TRIMPOLYA 1

#define BASEA 1
#define BASEC 2
#define BASEG 4
#define BASET 8

unsigned char *base_tab_init()
{
    unsigned char *t = malloc(sizeof(char)*256);
    memset(t, 0, sizeof(char)*256);
    t['A'] = BASEA;
    t['C'] = BASEC;
    t['G'] = BASEG;
    t['T'] = BASET; 
    t['a'] = BASEA;
    t['c'] = BASEC;
    t['g'] = BASEG;
    t['t'] = BASET;
    return t;
}
unsigned char *base_tab_rev_init()
{
    unsigned char *t = malloc(sizeof(char)*256);
    memset(t, 0, sizeof(char)*256);
    t['A'] = BASET;
    t['C'] = BASEG;
    t['G'] = BASEC;
    t['T'] = BASEA; 
    t['a'] = BASET;
    t['c'] = BASEG;
    t['g'] = BASEC;
    t['t'] = BASEA;
    return t;
}

struct trimstat {
    uint64_t all_fragments;
    uint64_t trimmed;
    uint64_t small;
    uint64_t dropped;
};

// Transposase recognition sequences
// 19bp Mosaic Ends: CTGTCTCTTATACACATCT
const char *me = "CTGTCTCTTATACACATCTGACGTC";
//                    CTGTCTCTTATACACATCTGACGTC
const char *rev_me = "TGTGTATAAGAGACAG";
const char *polyA = "AAAAAAAAAAAAAAAAA";
const char *polyT = "TTTTTTTTTTTTTTTTT";

// AGCGTCAGATGTGTATAAGAGACAG
const char *code2seq = "NACNGNNNTNN";

#define MINI_SKIP 1
#define TRIMMED   2
#define DROP_POLL 3

int countbits(uint64_t x)
{
    int i;
    int l = 0;
    for (i = 0; i < 16; i++) {
        if (((x>>(i*4))&0xf) > 0) l++;
    }
    return l;
}
static int usage()
{
    fprintf(stderr, "* This program used to trim the mosic ends and adaptor sequence introduced at Tn5 library.\n");
    fprintf(stderr, "Usage: Dyncut [options] read1.fq.gz [read2.fq.gz]\n");
    fprintf(stderr,"\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr,"  -1 [output1.fq]      Output fastq 1 file or smart pairing file (ignore fastq 2).\n");
    fprintf(stderr,"  -2 [output2.fq]      Output fastq 2 file.\n");
    fprintf(stderr,"  -t [1]               Threads.\n");
    fprintf(stderr,"  -r [10000000]        Records cached per chunk.\n");
    fprintf(stderr,"  -p                   Smart pairing.\n");
    fprintf(stderr,"  -mode [Tn5|polyA]    Trim mode, Tn5 for mosic ends, ployA for ployA tails.\n");
    fprintf(stderr,"  -l [20]              Minimal fragment to keep.\n");
    fprintf(stderr,"  -m [1]               Allowed mismatches.\n");
    fprintf(stderr,"  -report [rept.json]  Export report summary.\n"); 
    fprintf(stderr,"  -order               Export records in order (According to input records).\n");
    fprintf(stderr,"  -d                   Drop if reversed ME sequence detected.\n");
    fprintf(stderr,"\n"
            "Notes:\n"
            "The program is designed to trim TN5 transposase introduced mosaic ends and\n"
            "adaptor populations. But it should also be work for other kinds of transposases\n"
            "by define a specific adaptor sequence.\n"
            "Demo 1:\n"
            " DyncutTn5 -tail 5 reads1.fq.gz reads2.fq.gz | bwa mem -Cpt8 ref.fa - \\\n"
            " | samtools view -Sb -|sambamba sort /dev/stdin -|samtools rmdup /dev/stdin aln.bam\n"
            "Demo 2:\n"
            " DyncutTn5 -tail 5 reads1.fq.gz reads2.fq.gz -1 out1.fq.gz -2 out2.fq.gz &&\\\n"
            " bowtie2 -x ref.fa -1 out1.fq.gz -2 out2.fq.gz | samtools view -Sb - | ...\n"
        );
    return 1;
}

struct encode {
    int l;
    uint64_t x;
};
struct encode *revcode(struct encode *c)
{
    struct encode *r = malloc(sizeof(*r));
    uint8_t x[9];
    memset(x, 0, 9);
    x[1] = 8;
    x[2] = 4;
    x[4] = 2;
    x[8] = 1;
    r->l = c->l;
    r->x = 0;
    int i;
    for (i = c->l-1; i >= 0; --i) r->x = (r->x<<4) | x[(c->x>>(i*4))&0xf];
    return r;
};
struct encode *str_encode(const char *s, unsigned char *tab)
{
    int i;
    int l = strlen(s);
    struct encode *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    // only encode first 16 bases or shorter
    x->l = l > 16 ? 16 : l;
    for (i = 0; i < x->l; ++i) {
        x->x = x->x<<4 | (tab[s[i]]&0xf);
    }
    return x;
}

static struct args {
    const char *input_fname1;
    const char *input_fname2;
    const char *output_fname1;
    const char *output_fname2;
    const char *report_fname;
    const char *fail_fname;
    const char *se_fname;

    int smart_pair;
    int n_thread;
    
    int skip_flag;
    int se_mode;
    int mismatch;
    int tail_3;
    int mini_frag;
    int rev_trimmed;    
    const char *adaptor;
    unsigned char *base_tab;
    unsigned char *rev_tab;

    FILE *r1_out; // change gzipped output to normal file, saving time..
    FILE *r2_out;
    
    FILE *report_fp;
    FILE *fail_fp;

    int is_pe;
    int chunk_size;

    int mode;    
    struct encode *me_or_ada;
    struct encode *revada;

    struct trimstat stat;

    struct fastq_handler *fastq;
    
} args = {
    .input_fname1 = NULL,
    .input_fname2 = NULL,
    .output_fname1 = NULL,
    .output_fname2 = NULL,
    .report_fname = NULL,
    .fail_fname = NULL,
    .se_fname = NULL,
    
    .smart_pair = 0,
    .n_thread = 1,
    //.no_parse_flag = 0,
    
    .skip_flag = 0,
    .se_mode = 0,
    .mismatch = 1,
    .mini_frag = 18,
    .tail_3 = 0,
    .rev_trimmed = 0,
    .adaptor = NULL,
    .base_tab = NULL,
    .report_fp = NULL,
    .fail_fp = NULL,

    .r1_out = NULL,
    .r2_out = NULL,

    .mode = TRIMME,
    .me_or_ada = NULL,
    .revada =NULL,
    .chunk_size = 10000000, // 10M
    .is_pe = 0,
    .stat = {0,0,0},

    .fastq = NULL,
};
static void memory_release()
{
    if (args.report_fp) {
        fprintf(args.report_fp, "{\n\t\"All fragments\":\"%llu\",\n", args.stat.all_fragments);
        fprintf(args.report_fp, "\t\"Trimmed fragments\":\"%llu\",\n", args.stat.trimmed);
        fprintf(args.report_fp, "\t\"Fragments smaller than %d\":\"%llu\",\n", args.mini_frag, args.stat.small);
        fprintf(args.report_fp, "\t\"Rev adaptor polluated fragments\":\"%llu\",\n", args.stat.dropped);
        uint64_t n_pass = args.stat.all_fragments - args.stat.small - args.stat.dropped;        
        fprintf(args.report_fp, "\t\"Reads passed QC (rate)\":\"%llu (%.2f%%)|\"\n", n_pass*(args.is_pe+1), (float)n_pass/args.stat.all_fragments*100);
        fprintf(args.report_fp, "}\n");
        fclose(args.report_fp);
    }
                
    free(args.base_tab);
    free(args.rev_tab);
    free(args.me_or_ada);
    free(args.revada);

    if (args.r1_out) fclose(args.r1_out);
    if (args.r2_out) fclose(args.r2_out);
}
static int parse_args(int argc, char **argv)
{
    int i;
    const char *mis = 0;
    const char *t3 = 0;
    const char *mini = 0;
    const char *thread = 0;
    const char *record = 0;
    const char *mode = 0;
    
    for (i = 1; i < argc;) {

        const char *a = argv[i++];
        const char **var = 0;
        
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return usage();
        //if (strcmp(a, "-v") == 0 || strcmp(a, "--version") == 0) return get_version();
        if (strcmp(a, "-1") == 0) var = &args.output_fname1;
        else if (strcmp(a, "-2") == 0) var = &args.output_fname2;
        else if (strcmp(a, "-adaptor") == 0) var = &args.adaptor;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-fail") == 0) var = &args.fail_fname;
        else if (strcmp(a, "-m") == 0) var = &mis;
        else if (strcmp(a, "-l") == 0) var = &mini;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tail") == 0) var = &t3;
        else if (strcmp(a, "-r") == 0) var = &record;
        else if (strcmp(a, "-mode") == 0) var = &mode;
        else if (strcmp(a, "-d") == 0) {
            args.rev_trimmed = 1;
            continue;
        }
        else if (strcmp(a, "-p") == 0) {
            args.smart_pair = 1;
            continue;
        }

        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1]) error("Unknown parameter. %s", a);
        
        if (args.input_fname1 == 0) {
            args.input_fname1 = a;
            continue;
        }

        if (args.input_fname2 == 0) {
            args.input_fname2 = a;
            continue;
        }

        error("Unknown argument: %s, use -h see help information.", a);
    }

    if (mis) args.mismatch = str2int((char*)mis);
    if (mini) args.mini_frag = str2int((char*)mini);
    if (thread) args.n_thread = str2int((char*)thread);
    if (t3) args.tail_3 = str2int((char*)t3);
    if (record) args.chunk_size = str2int((char*)record);
    if (mode) {
        if (strcmp(mode, "Tn5")==0) args.mode = TRIMME;
        else if (strcmp(mode, "polyA") == 0) args.mode = TRIMPOLYA;
        else error("Unknown mode, %s", mode);
    }
    
    args.base_tab = base_tab_init();
    args.rev_tab = base_tab_rev_init();

    if (args.mode == TRIMME) {
        args.me_or_ada = str_encode(me, args.base_tab);
        args.revada = str_encode(rev_me, args.base_tab);
    }
    else if (args.mode == TRIMPOLYA) {
        args.me_or_ada = str_encode(polyA, args.base_tab);
        args.revada = str_encode(polyT, args.base_tab);    
    }
             
    if (args.input_fname1 == NULL && (!isatty(fileno(stdin))))
        args.input_fname1 = "-";    
    if (args.input_fname1 == NULL) error("Fastq file(s) must be set!");

    args.fastq = fastq_handler_init(args.input_fname1, args.input_fname2, args.smart_pair, args.chunk_size);

    int state = fastq_handler_state(args.fastq);
    if (state == FH_SMART_PAIR || state == FH_PE ) args.is_pe = 1;

    if (args.output_fname1 != NULL) {
        args.r1_out = fopen(args.output_fname1, "w");
        if (args.r1_out == NULL) error("%s : %s.", args.output_fname1, strerror(errno));
        if (args.output_fname2 != NULL) {
            args.r2_out = fopen(args.output_fname2, "w");
            if (args.r2_out == NULL) error("%s : %s.", args.output_fname2, strerror(errno));
        }
    }

    if (args.report_fname) {
        args.report_fp = fopen(args.report_fname, "w");
        if (args.report_fp == NULL) error("%s : %s.", args.report_fname, strerror(errno));
    }

    if (args.fail_fname) {
        args.fail_fp = fopen(args.fail_fname, "w");
        if (args.fail_fp == NULL) error("%s : %s.", args.fail_fname, strerror(errno));
    }

    return 0;
}
int find_sequence_adaptor(const char *s, const struct encode *a, int m, unsigned char const *tab)
{
    int l, n = 0;
    int i;
    l = strlen(s);    
    uint64_t x = 0;
    //uint64_t mask = 0xFFFFFFFFFFFFFFFF;
    for (i = 0; i < l; ++i) {
        x = x<<4|(tab[s[i]]&0xf);
        if (++n >= a->l) {
            if (x == a->x ||countbits(x&a->x) >= a->l - m) return i-a->l+1;
        }
    }
    // tails
    int mis = m;
    uint64_t x1 = a->x;
    for (i = 1; i < a->l-3; ++i) {
        mis = mis <= 0 ? 0 : --mis;
        int l1 = a->l - i;
        x1 = x1>>4;
        if (countbits(x&x1)>= l1 -mis) return l - a->l + i+1;
    }
    return -1;
}
int find_sequence_adaptor_rev(const char *s, int l, const struct encode *a, int m, unsigned char const *tab)
{
    int n = 0;
    int i;
    uint64_t x = 0;
    for (i = 0; i < l; ++i) {
        x = x<<4|(tab[s[i]]&0xf);
        if (++n >= a->l) {
            if (x== a->x || countbits(x&a->x) >= a->l -m) return i+1;
        }
    }

    return -1;
}
void trim_3end(struct args *opts, struct bseq_pool *p)
{
    int i;
    if (opts->is_pe) {
        for (i = 0; i < p->n; ++i) {
            struct bseq *b = &p->s[i];
            int l0, l1;
                            
            l0 = find_sequence_adaptor(b->s0, opts->me_or_ada, opts->mismatch, opts->base_tab);
            l1 = find_sequence_adaptor(b->s1, opts->me_or_ada, opts->mismatch, opts->base_tab);

            if (b->l0 - l0 < 5 && l1 == -1) l0 = -1;
            if (b->l1 - l1 < 5 && l0 == -1) l1 = -1;
            // no ME fragment found at both ends
            if (l1 == -1 && l0 == -1) {
                if (opts->tail_3 > 0 && b->l1 > opts->tail_3) {
                    b->l1 -= opts->tail_3;
                    b->l0 -= opts->tail_3;
                    if (b->l1 < args.mini_frag) b->flag = MINI_SKIP;
                }
                // continue;
            }
            
            // consider PE reads are not same length, treat as SEs??

            // ME found
            else if (l0 == l1) {
                b->l0 = l0;
                b->l1 = l1;
                if (b->l1 < args.mini_frag) b->flag = MINI_SKIP;
                else b->flag = TRIMMED;
                // continue;
            }

            // Inconsistant, trim as many as possible
            else if (l0 != l1) {
                // equal length
                if (l0 > 0 && l1 > 0) {
                    b->l0 = l0 > l1 && l1 > 0 ? l1 : l0;
                    b->l1 = l1 > l0 && l0 > 0 ? l0 : l1;
                    if (b->l1 < args.mini_frag) b->flag = MINI_SKIP;
                    else b->flag = TRIMMED;
                }
            }
            b->s0[b->l0] = 0;
            b->s1[b->l1] = 0;
            if (b->q0) b->q0[b->l0] = 0;
            if (b->q1) b->q1[b->l1] = 0;
        }
    }
    else { // Single end
        for (i = 0; i < p->n; ++i) {
            struct bseq *b = &p->s[i];
            int l0;
            l0 = find_sequence_adaptor(b->s0, opts->me_or_ada, opts->mismatch, opts->base_tab);

            if (b->l0 - l0 < 5) l0 = -1;
            // no ME fragment found at both ends
            if (l0 == -1) {
                if (opts->tail_3 > 0 && b->l0 > opts->tail_3) {
                    b->l0 -= opts->tail_3;
                    if (b->l0 < args.mini_frag) b->flag = MINI_SKIP;
                }
                // continue;
            }
            else {
                b->l0 = l0;
                if (b->l0 < args.mini_frag) b->flag = MINI_SKIP;
                else b->flag = TRIMMED;
            }
            b->s0[b->l0] = 0;
            if (b->q0) b->q0[b->l0] = 0;
        }
    }

}
void trim_5end(struct args *opts, struct bseq_pool *p)
{
    int l0, l1;
    int i;
    if (opts->rev_trimmed) { // drop all reversed adaptor polluation instead of trim
        for (i = 0; i < p->n; ++i) {
            struct bseq *b = &p->s[i];

            l0 = find_sequence_adaptor_rev(b->s0, b->l0, opts->revada, opts->mismatch, opts->base_tab);
            if (l0 > 0) {
                //b->l0 -= l0;
                //memmove(b->s0, b->s0 + l0, b->l0);
                b->flag = DROP_POLL;
            }
            else if (b->l1 > 0) {
                l1 = find_sequence_adaptor_rev(b->s1, b->l1, opts->revada, opts->mismatch, opts->base_tab);
                if (l1 > 0) {
                    //b->l1 -= l1;
                    //memmove(b->s1, b->s1 + l1, b->l1);
                    b->flag = DROP_POLL;
                }
            }

        }
    }
}
static void write_out(void *_data)
{
    struct bseq_pool *p = (struct bseq_pool*)_data;
    struct args *opts = p->opts;

    FILE *fp1 = opts->r1_out == NULL ? stdout : opts->r1_out;
    FILE *fp2 = opts->r2_out == NULL ? fp1 : opts->r2_out;
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        opts->stat.all_fragments++;
        
        if (b->flag == MINI_SKIP) {
            opts->stat.small++;
        }
        else if (b->flag == DROP_POLL) {
            opts->stat.dropped++;
        }
        else {
            if (b->flag == TRIMMED) opts->stat.trimmed++;
            fprintf(fp1, "%c%s\n%s\n", b->q0 ? '@' : '>', b->n0, b->s0);
            if (b->q0) fprintf(fp1, "+\n%s\n", b->q0);
            if (b->l1 > 0) {
                fprintf(fp2, "%c%s\n%s\n", b->q1 ? '@' : '>', b->n0, b->s1);
                if (b->q1) fprintf(fp2, "+\n%s\n", b->q1);
            }
        }        
    }   
    bseq_pool_destroy(p);
}

static void *run_it(void *_p, int idx)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    struct args *opts = p->opts;
    trim_3end(opts, p);
    trim_5end(opts, p);
    return p;
}

int fastq_trim_adaptors(int argc, char **argv)
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
            block = thread_pool_dispatch2(p, q, run_it, b, 0);
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
