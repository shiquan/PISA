#include "utils.h"
#include "fastq.h"
#include "dict.h"
#include "htslib/thread_pool.h"
#include "htslib/kstring.h"
#include "number.h"
#include "read_tags.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static struct args {
    const char *r1_fname;
    const char *r2_fname;
    struct fastq_handler *fastq;
    const char *report_fname;
    FILE *fp_report;

    const char *out1_fname;
    const char *out2_fname;
    FILE *fp_out1;
    FILE *fp_out2;
    
    const char *parse_rules;
    int n_bc;
    struct bc_reg *bcs;

    struct bc_reg *r1;
    struct bc_reg *r2;

    int n_thread;
    int smart_pair;
    int qual_thres;
    int dropN;
    int chunk_size;

    int no_warnings;
    // stats of reads
    uint64_t raw_reads;
    uint64_t reads_pass_qc;
    uint64_t barcode_exactly_matched;
    uint64_t filtered_by_barcode;
    uint64_t filtered_by_lowqual;
} args = {
    .r1_fname    = NULL,
    .r2_fname    = NULL,
    .fastq       = NULL,
    .report_fname= NULL,
    .fp_report   = NULL,
    .out1_fname  = NULL,
    .out2_fname  = NULL,
    .fp_out1     = NULL,
    .fp_out2     = NULL,
    
    .parse_rules = NULL,
    .n_bc        = 0,
    .bcs         = NULL,
    .r1          = NULL,
    .r2          = NULL,
    .smart_pair  = 0,
    .qual_thres  = 20,
    .dropN       = 0,
    .chunk_size  = 1000000,
    .n_thread    = 4,
    .raw_reads   = 0,
    .reads_pass_qc = 0,
    .barcode_exactly_matched = 0,
    .filtered_by_barcode = 0,
    .filtered_by_lowqual = 0,
};

struct bc_reg {
    int rd;
    int st; // 0 offset
    int ed; // 1 offset
    char *raw_tag;
    char *corr_tag;
    struct dict *wl;
};

struct dict *read_wl(const char *fn, int mis)
{
    struct dict *wl = dict_init();
    if ( dict_read(wl, fn) ) {
        dict_destroy(wl);
        error("Failed to read barcodes from %s.", fn);
    }

    int i;
    kstring_t str = {0,0,0};
    int n = dict_size(wl); // current keys of white list
    for (i = 0; i < n; ++i) {
        str.l = 0;
        kputs(dict_name(wl, i), &str);

        if (mis == 0) continue;
        
        int j;
        for (j = 0; j < str.l; ++j) {
            char o = str.s[j];
            str.s[j] = 'N';
            int idx = dict_query2(wl, str.s);
            if (idx == -2) continue; // already deleted
            if (idx >= 0) {
                if (args.no_warnings == 0)
                    warnings("Barcode %s has multi hits in the white list.", str.s);

                dict_del(wl, str.s);
                str.s[j] = o;
                continue;
            }
            dict_push2(wl, str.s, i);
            str.s[j] = o;
        }
    }
    free(str.s);
    return wl;
}

char *bseq_subset_seq(struct bseq *b, struct bc_reg *r)
{
    if (r == NULL) return NULL;
    kstring_t str = {0,0,0};
    if (r->rd == 1) {
        if (r->st == -1) kputs(b->s0.s, &str);
        else {
            kputsn(b->s0.s+r->st, r->ed - r->st, &str);
            kputs("", &str);
        }
    } else if (r->rd == 2) {
        if (r->st == -1) kputs(b->s1.s, &str);
        else {
            kputsn(b->s1.s+r->st, r->ed - r->st, &str);
            kputs("", &str);
        }
    } else {
        error("Only support 2 reads now.");
    }
    return str.s;
}
char *bseq_subset_qual(struct bseq *b, struct bc_reg *r)
{
    if (r == NULL) return NULL;
    kstring_t str = {0,0,0};
    if (r->rd == 1) {
        if (b->q0.l == 0 || b->q0.l < r->st) return NULL;        
        if (r->st == -1) {
            kputs(b->q0.s, &str);
        } else if (r->ed > r->st){
            kputsn(b->q0.s+r->st, r->ed - r->st, &str);
            kputs("", &str);
        } else {
            kputs(b->q0.s+r->st, &str);
        }
    } else if (r->rd == 2) {
        if (b->q1.l == 0) return NULL;
        if (r->st == -1) {
            kputs(b->q1.s, &str);
        } else if (r->ed > r->st) {
            kputsn(b->q1.s+r->st, r->ed - r->st, &str);
            kputs("", &str);
        } else {
            kputs(b->q1.s+r->st, &str);
        }
    } else {
        error("Only support 2 reads now.");
    }
    return str.s;
}
char *correct_bc(struct bc_reg *r, const char *val, int *exact)
{
    *exact = 0;
    if (r->wl == NULL) return NULL;
    int idx = dict_query(r->wl, val);
    if (idx >= 0) {
        *exact = 1;
        return dict_name(r->wl, idx); // Exactly match
    }

    int i;
    int l = strlen(val);
    kstring_t str= {0, 0, 0};
    kputs(val, &str);
    for (i = 0; i < l; ++i) {
        str.s[i] = 'N';
        idx = dict_query(r->wl, str.s);
        if (idx > 0) {
            free(str.s);
            return dict_name(r->wl, idx);
        }
    }
    free(str.s);
    return NULL;
}
int parse_region(const char *_s, struct bc_reg *r)
{
    r->st = -1;
    r->ed = -1;
    
    char *s0 = strdup(_s);
    char *s = s0;
    if (s[0] != 'R') { free(s0); return 1;}
    if (s[1] == '1') r->rd = 1;
    else if (s[1] == '2') r->rd = 2;
    else { free(s0); return 1; }

    if (s[2] == '\0') { free(s0); return 0; }
    s = s+2;

    if (*s == ',' || *s == ';') { free(s0); return 0; }
    if (*s != ':') { free(s0); return 1; }

    s = s+1;

    char *p = s;
    int n = 0;
    for (; s && isdigit(*s); s++) n++;
    if (*s == '-') {
        s++;
        p[n] = '\0';   
    }
    r->st = str2int(p) - 1;

    if (s) r->ed = str2int(s);
    free(s0);
    return 0;
}
// CB,R1:1-10,whitelist.txt,CB,1;R1,R1:11-60;R2,R2
void parse_rules(const char *rule)
{
    kstring_t str = {0,0,0};

    kputs(rule, &str);
    int n;
    int *s = ksplit(&str, ';', &n);

    int alloc = n -1;
    if (alloc == 0)  error("Unrecognised rule format.");

    args.bcs = malloc(sizeof(struct bc_reg)*alloc);
    int n0 = 0;
    
    int i;
    for (i = 0; i < n; ++i) {
        char *ru = str.s+s[i];
        if (strlen(ru) < 2)  error("Unrecognised rule format.");
        
        if (ru[0] == 'R' && ru[1] == '1') {
            if(ru[2] != ',') error("Unrecognised rule format.");
            ru = ru + 3;
            args.r1 = malloc(sizeof(struct bc_reg));
            memset(args.r1, 0, sizeof(struct bc_reg));
            if (parse_region(ru, args.r1)) error("Unrecognised rule format.");
        } else if (ru[0] == 'R' && ru[1] == '2') {
            if(ru[2] != ',') error("Unrecognised rule format.");
            ru = ru + 3;
            args.r2 = malloc(sizeof(struct bc_reg));
            memset(args.r2, 0, sizeof(struct bc_reg));
            if (parse_region(ru, args.r2)) error("Unrecognised rule format.");
        } else {
            if (ru[2] != ',') error("Unrecognised rule format.");
            struct bc_reg *r0 = &args.bcs[n0++];
            memset(r0, 0, sizeof(struct bc_reg));
            kstring_t temp = {0,0,0};
            kputs(ru, &temp);
            int n1;
            int *s0 = ksplit(&temp, ',', &n1);
            r0->raw_tag = strdup(temp.s);
            if (parse_region(temp.s+s0[1], r0)) error("Unrecognised rule format.");

            int mis = 0;
            if (n1 >=4) {
                if (strlen(temp.s+s0[3]) != 2) error("Unrecognised rule format.");
                if (*(temp.s+s0[4]) == '1') mis = 1;
                r0->wl = read_wl(temp.s+s0[2], mis);
                r0->corr_tag = strdup(temp.s+s0[3]);
            }

            free(temp.s);
            free(s0);
        }
    }

    args.n_bc = n0;
    free(s);
    free(str.s);
}
static int parse_args(int argc, char **argv)
{
    if ( argc == 1 ) return 1;
    
    int i;
    const char *thread = NULL;
    const char *qual_thres = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        
        if (strcmp(a, "-1") == 0) var = &args.out1_fname;
        else if (strcmp(a, "-2") == 0) var = &args.out2_fname;
        else if (strcmp(a, "-rule") == 0) var = &args.parse_rules;
        else if (strcmp(a, "-t") == 0) var = &thread; 
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;       
        else if (strcmp(a, "-p") == 0) {
            args.smart_pair = 1;
            continue;
        }
        else if (strcmp(a, "-dropN") == 0) {
            args.dropN = 1;
            continue;
        }
        else if (strcmp(a, "-nw") == 0) {
            args.no_warnings = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (a[0] == '-' && a[1]) error("Unknown parameter %s.", a);
        if (args.r1_fname == 0) {
            args.r1_fname = a;
            continue;
        }
        if (args.r2_fname == 0) {
            args.r2_fname = a;
            continue;
        }
        error("Unknown argument: %s, use -h see help information.", a);
    }

    if (args.parse_rules == NULL) error("Option -rule is required.");

    parse_rules(args.parse_rules);
    
    if (thread) args.n_thread = str2int((char*)thread);
    if (qual_thres) {
        args.qual_thres = str2int((char*)qual_thres);
        LOG_print("Average quality below %d will be drop.", args.qual_thres);
    }

    if (args.r1_fname == NULL && (!isatty(fileno(stdin)))) args.r1_fname = "-";
    if (args.r1_fname == NULL) error("Fastq file(s) must be set.");
        
    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        CHECK_EMPTY(args.fp_report, "%s : %s.", args.report_fname, strerror(errno));
    }
    else args.fp_report = stderr;

    if (args.out1_fname) {
        args.fp_out1 = fopen(args.out1_fname, "w");
        if (args.fp_out1 == NULL) error("%s: %s.", args.out1_fname, strerror(errno));
        if (args.out2_fname) {
            args.fp_out2 = fopen(args.out2_fname,"w");
            if (args.fp_out2 == NULL) error("%s: %s.", args.out2_fname, strerror(errno));
        }
    } else {
        args.fp_out1 = stdout;
    }
    
    args.fastq = fastq_handler_init(args.r1_fname, args.r2_fname, args.smart_pair, args.chunk_size);
    if (args.fastq == NULL) error("Failed to init input fastq.");
    
    return 0;
}

static void memory_release()
{
    if (args.fp_out1 != stdout) fclose(args.fp_out1);
    if (args.fp_out2) fclose(args.fp_out2);
    if (args.fp_report && args.fp_report != stderr) fclose(args.fp_report);
    int i;
    for (i = 0; i < args.n_bc; ++i) {
        struct bc_reg *r = &args.bcs[i];
        if (r->wl) dict_destroy(r->wl);
        if (r->raw_tag) free(r->raw_tag);
        if (r->corr_tag) free(r->corr_tag);
    }
    
    free(args.bcs);

    struct bc_reg *r = args.r1;
    if (r->wl) dict_destroy(r->wl);
    if (r->raw_tag) free(r->raw_tag);
    if (r->corr_tag) free(r->corr_tag);
    free(args.r1);
    
    if (args.r2) {
        struct bc_reg *r = args.r2;
        if (r->wl) dict_destroy(r->wl);
        if (r->raw_tag) free(r->raw_tag);
        if (r->corr_tag) free(r->corr_tag);
        free(args.r2);        
    }
    fastq_handler_destory(args.fastq);
}
static int write_report()
{
    if (args.fp_report == NULL) return 1;
    fprintf(args.fp_report, "Number of Fragments,%"PRIu64"\n", args.raw_reads);
    fprintf(args.fp_report, "Fragments pass QC,%"PRIu64"\n", args.reads_pass_qc);
    fprintf(args.fp_report, "Fragments with Exactly Matched Barcodes,%"PRIu64"\n", args.barcode_exactly_matched);
    fprintf(args.fp_report, "Fragments with Failed Barcodes,%"PRIu64"\n", args.filtered_by_barcode);
    // fprintf(args.fp_report, "Fragments Filtered on Low Quality,%"PRIu64"\n", args.filtered_by_lowqual);
    // fprintf(args.fp_report, "Q30 bases in Reads,%.1f%%\n", (float)args.q30_bases_reads/(args.bases_reads+1)*100);
    
    if (args.fp_report != stderr) fclose(args.fp_report);
    return 0;
}
static void *run_it(void *_p)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        int j;
        for (j = 0; j < args.n_bc; ++j) {
            struct bc_reg *r = &args.bcs[j];
            char *val = bseq_subset_seq(b, r);
            char *name1 = fname_update_tag(b->n0.s, r->raw_tag, val);
            b->n0.l = 0;
            kputs(name1, &b->n0);
            free(name1);
            if (r->corr_tag) {
                int ex;
                char *val0 = correct_bc(r, val, &ex);
                if (ex) b->flag = FQ_FLAG_BC_EXACTMATCH;
                if (val0 == NULL) {
                    b->flag = FQ_FLAG_BC_FAILURE;
                } else {
                    name1 = fname_update_tag(b->n0.s,r->corr_tag, val0);
                    b->n0.l = 0;
                    kputs(name1, &b->n0);
                    free(name1);
                }
            }
            free(val);
        }

        char *r1 = NULL;
        char *q1 = NULL;
        char *r2 = NULL;
        char *q2 = NULL;
        r1 = bseq_subset_seq(b, args.r1);
        q1 = bseq_subset_qual(b, args.r1);
        if (r1 == NULL) error("Empty read one.");
        r2 = bseq_subset_seq(b, args.r2);
        q2 = bseq_subset_qual(b, args.r2);

        // update reads
        b->s0.l = 0;
        b->q0.l = 0;
        kputs(r1, &b->s0);
        if (q1) kputs(q1, &b->q0);

        b->s1.l = 0;
        b->q1.l = 0;
        if (r2) kputs(r2, &b->s1);
        if (q2) kputs(q2, &b->q1);

        if (r1) free(r1);
        if (q1) free(q1);
        if (r2) free(r2);
        if (q2) free(q2);
    }
    return p;
}

static void write_out(void *_p)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    
    FILE *fp1 = args.fp_out1 == NULL ? stdout : args.fp_out1;
    FILE *fp2 = args.fp_out2 == NULL ? fp1 : args.fp_out2;

    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        
        args.raw_reads++;
        
        if (b->flag == FQ_FLAG_BC_FAILURE) {
            args.filtered_by_barcode++;
            continue;
        }

        if (b->flag == FQ_FLAG_READ_QUAL) {
            args.filtered_by_lowqual++;
            continue; // just skip ALL low quality reads
        }

        if (b->flag == FQ_FLAG_BC_EXACTMATCH) {
            args.barcode_exactly_matched++;
        }

        args.reads_pass_qc++;
        fprintf(fp1, "%c%s\n%s\n", b->q0.l ? '@' : '>', b->n0.s, b->s0.s);
        if (b->q0.l) fprintf(fp1, "+\n%s\n", b->q0.s);
        if (b->s1.l > 0) {
            fprintf(fp2, "%c%s\n%s\n", b->q1.l ? '@' : '>', b->n0.s, b->s1.s);
            if (b->q1.l) fprintf(fp2, "+\n%s\n", b->q1.s);
        }
    }
    bseq_pool_destroy(p);
    fflush(fp1);
    if (fp2 != fp1) fflush(fp2);
}

extern int fastq_parse2_usage();

int fastq_parse3(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return fastq_parse2_usage();

        
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        struct bseq_pool *pool = fastq_read(args.fastq, NULL);
        
        if (pool == NULL) break;
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, pool, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct bseq_pool *pool = (struct bseq_pool*)hts_tpool_result_data(r);
                write_out(pool);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);
    
    while ((r = hts_tpool_next_result(q))) {
        struct bseq_pool *pool = (struct bseq_pool*)hts_tpool_result_data(r);
        write_out(pool);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    write_report();
    memory_release();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);

    return 0;    
}


int fastq_parse2(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return fastq_parse2_usage();

    
    struct bseq_pool *b;

#pragma omp parallel private(b) num_threads(args.n_thread)
    for (;;) {
        
#pragma omp critical (read)
        b = fastq_read(args.fastq, NULL);
        if (b == NULL) break;
        b = run_it(b);
        
#pragma omp critical (write)
        write_out(b);            
    }
    
    write_report();
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}