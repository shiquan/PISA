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
    const char *r3_fname;
    const char *r4_fname;
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
    int order;
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
    .r3_fname    = NULL,
    .r4_fname    = NULL,
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
    .order       = 0,
    .dropN       = 0,
    .chunk_size  = 1000,
    .n_thread    = 4,
    .raw_reads   = 0,
    .reads_pass_qc = 0,
    .barcode_exactly_matched = 0,
    .filtered_by_barcode = 0,
    .filtered_by_lowqual = 0,
};

struct bc_reg0 {
    int rd;
    int st; // 0 offset
    int ed; // 1 offset
    struct dict *wl;
};

struct bc_reg {
    int n;
    char *raw_tag;
    char *corr_tag;
    struct bc_reg0 *r;
};
struct dict *build_mis(struct dict *wl)
{
    int i;
    kstring_t str = {0,0,0};
    int n = dict_size(wl); // current keys of white list
    for (i = 0; i < n; ++i) {
        str.l = 0;
        kputs(dict_name(wl, i), &str);
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
struct dict *read_wl(const char *fn, int mis)
{
    struct dict *wl = dict_init();
    if ( dict_read(wl, fn, 0) ) {
        dict_destroy(wl);
        error("Failed to read barcodes from %s.", fn);
    }

    if (mis == 0) return wl;

    return build_mis(wl);
}
struct dict *read_wl_cached(const char **bcs, int l, int mis)
{
    struct dict *wl = dict_init();
    int i;
    for (i = 0; i < l;i++) {
        // if (bcs[i] == NULL) continue;
        dict_push(wl, bcs[i]);
    }
    if (mis == 0) return wl;
    return build_mis(wl);
}

void bseq_subset_str0(kstring_t *from, kstring_t *to, int st, int ed)
{
    if (st == -1) kputs(from->s, to);
    else if (ed == -1) kputs(from->s+st, to);
    else {
        kputsn(from->s+st, ed - st, to);
        kputs("", to);
    }
}
char *bseq_subset_seq(struct bseq *b, int rd, int st, int ed)
{
    kstring_t str = {0,0,0};
    
    if (rd == 1) {
        kstring_t *from = &b->s0;
        bseq_subset_str0(from, &str, st, ed);
    } else if (rd == 2) {
        kstring_t *from = &b->s1;
        bseq_subset_str0(from, &str, st, ed);
    } else if (rd == 3) {
        kstring_t *from = &b->s2;
        bseq_subset_str0(from, &str, st, ed);
    } else if (rd == 4) {
        kstring_t *from = &b->s3;
        bseq_subset_str0(from, &str, st, ed);        
    } else {
        error("Only support 1-4 read files now.");
    }
    return str.s;
}
char *bseq_subset_qual(struct bseq *b, int rd, int st, int ed)
{
    kstring_t str = {0,0,0};
    if (rd == 1) {
        if (b->q0.l == 0) return NULL; 
        kstring_t *from = &b->q0;
        bseq_subset_str0(from, &str, st, ed);
    } else if (rd == 2) {
        if (b->q1.l == 0) return NULL; 
        kstring_t *from = &b->q1;
        bseq_subset_str0(from, &str, st, ed);
    } else if (rd == 3) {
        if (b->q2.l == 0) return NULL; 
        kstring_t *from = &b->q2;
        bseq_subset_str0(from, &str, st, ed);
    } else if (rd == 4) {
        if (b->q3.l == 0) return NULL; 
        kstring_t *from = &b->q3;
        bseq_subset_str0(from, &str, st, ed);
    } else {
        error("Only support 1-4 read files now.");
    }
    return str.s;
}
char *correct_bc(struct dict *wl, const char *val, int *exact)
{
    *exact = 0;
    int idx = dict_query(wl, val);
    if (idx >= 0) {
        *exact = 1;
        return dict_name(wl, idx); // Exactly match
    } else {
        *exact = 0;
        int j;
        int l = strlen(val);
        kstring_t str= {0, 0, 0};
        kputs(val, &str);
        for (j = 0; j < l; ++j) {
            char o = str.s[j];
            str.s[j] = 'N';
            idx = dict_query2(wl, str.s);
            if (idx > 0) {
                free(str.s);
                return dict_name(wl, idx);
            }
            str.s[j] = o;
        }
        free(str.s);        
    }
    return NULL;
}
int parse_region(const char *_s, struct bc_reg *r)
{
    struct bc_reg0 *r0 = &r->r[0];
    
    r0->st = -1;
    r0->ed = -1;
    
    char *s0 = strdup(_s);
    char *s = s0;
    if (s[0] != 'R') { free(s0); return 1;}
    if (s[1] == '1') r0->rd = 1;
    else if (s[1] == '2') r0->rd = 2;
    else if (s[1] == '3') r0->rd = 3;
    else if (s[1] == '4') r0->rd = 4;
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
    r0->st = str2int(p) - 1;

    if (s) r0->ed = str2int(s);
    free(s0);
    return 0;
}
// return 0 on add new bcs, return 1 on merged
int merge_bcs(struct bc_reg *bcs, int n0)
{
    if (n0 == 0) return 0;
    struct bc_reg *b0 = &bcs[n0];
    int i;
    for (i = 0; i < n0; ++i) {
        struct bc_reg *b = &bcs[i];       
        if (strcmp(b->raw_tag, b0->raw_tag) == 0) {
            b->r = realloc(b->r, (b->n+1)*sizeof(struct bc_reg0));
            b->r[b->n].rd = b0->r[0].rd;
            b->r[b->n].st = b0->r[0].st;
            b->r[b->n].ed = b0->r[0].ed;
            b->r[b->n].wl = b0->r[0].wl;
            b->n++;

            if (b0->raw_tag) free(b0->raw_tag);
            if (b0->corr_tag) free(b0->corr_tag);
            free(b0->r);
            return 1;
        }
    }
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
            args.r1->n = 1;
            args.r1->r = malloc(sizeof(struct bc_reg0));
            memset(args.r1->r, 0, sizeof(struct bc_reg0));
            if (parse_region(ru, args.r1)) error("Unrecognised rule format.");
        } else if (ru[0] == 'R' && ru[1] == '2') {
            if(ru[2] != ',') error("Unrecognised rule format.");
            ru = ru + 3;
            args.r2 = malloc(sizeof(struct bc_reg));
            memset(args.r2, 0, sizeof(struct bc_reg));
            args.r2->n = 1;
            args.r2->r = malloc(sizeof(struct bc_reg0));
            memset(args.r2->r, 0, sizeof(struct bc_reg0));
            if (parse_region(ru, args.r2)) error("Unrecognised rule format.");
        } else {
            if (ru[2] != ',') error("Unrecognised rule format.");
            struct bc_reg *r = &args.bcs[n0];
            memset(r, 0, sizeof(struct bc_reg));
            r->n = 1;
            r->r = malloc(sizeof(struct bc_reg0));
            memset(r->r, 0, sizeof(struct bc_reg0));
            kstring_t temp = {0,0,0};
            kputs(ru, &temp);
            int n1;
            int *s0 = ksplit(&temp, ',', &n1);
            r->raw_tag = strdup(temp.s);
            if (parse_region(temp.s+s0[1], r)) error("Unrecognised rule format.");

            int mis = 0;
            if (n1 >=4) {                
                if (strlen(temp.s+s0[3]) != 2) error("Unrecognised rule format.");
                if (*(temp.s+s0[4]) == '1') mis = 1;
                r->corr_tag = strdup(temp.s+s0[3]);
                struct bc_reg0 *r0 = &r->r[0];
                r0->wl = read_wl(temp.s+s0[2], mis);
            }

            free(temp.s);
            free(s0);
            if (merge_bcs(args.bcs, n0)==0) n0++; // n0 is real number of tags
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
    const char *code = NULL;
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
        else if (strcmp(a, "-x") == 0) var = &code;
        else if (strcmp(a, "-order") == 0) {
            args.order = 1;
            continue;
        }
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
        if (args.r3_fname == 0) {
            args.r3_fname = a;
            continue;
        }
        if (args.r4_fname == 0) {
            args.r4_fname = a;
            continue;
        }
        error("Unknown argument: %s, use -h see help information.", a);
    }

    if (args.parse_rules == NULL && code == NULL) error("Option -rule or -x is required.");
    if (args.parse_rules) {
        parse_rules(args.parse_rules);
    } else {
        if (strcmp(code, "C4") == 0) {

#include "DNBelabC4.h"            
            args.n_bc = 2;
            args.bcs = malloc(sizeof(struct bc_reg)*args.n_bc);
            struct bc_reg *cb = &args.bcs[0];
            memset(cb, 0, sizeof(*cb));
            
            cb->raw_tag = strdup("CR");
            cb->corr_tag = strdup("CB");
            cb->n = 2;
            cb->r = malloc(cb->n*sizeof(struct bc_reg0));
            struct bc_reg0 *s1 = &cb->r[0];
            struct bc_reg0 *s2 = &cb->r[1];
            s1->rd = s2->rd = 1;
            s1->st = 0;
            s1->ed = 10;
            s2->st = 10;
            s2->ed = 20;
            s1->wl = read_wl_cached(C4_barcodes, 1536, 1);
            s2->wl = read_wl_cached(C4_barcodes, 1536, 1);
            
            struct bc_reg *ub = &args.bcs[1];
            memset(ub, 0, sizeof(*ub));
            
            ub->n = 1;
            ub->raw_tag = strdup("UR");
            ub->r = malloc(sizeof(struct bc_reg0));
            ub->r->rd = 1;
            ub->r->st = 20;
            ub->r->ed = 30;
            ub->r->wl = NULL;
            
            args.r1 = malloc(sizeof(struct bc_reg));
            memset(args.r1, 0, sizeof(struct bc_reg));
            args.r1->n = 1;
            args.r1->r = malloc(sizeof(struct bc_reg0));
            struct bc_reg0 *r0 = &args.r1->r[0];
            memset(r0, 0, sizeof(*r0));
            
            r0->rd = 2;
            r0->st = 0;
            r0->ed = -1; // to the end                
        } else {
            error("Unknown code. %s", code);
        }
    }
    
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
    
    args.fastq = fastq_handler_init(args.r1_fname, args.r2_fname, args.r3_fname, args.r4_fname, args.smart_pair, args.chunk_size);
    if (args.fastq == NULL) error("Failed to init input fastq.");
    
    return 0;
}

static void memory_release()
{
    if (args.fp_out1 != stdout) fclose(args.fp_out1);
    if (args.fp_out2) fclose(args.fp_out2);

    int i;
    for (i = 0; i < args.n_bc; ++i) {
        struct bc_reg *r = &args.bcs[i];
        int j = 0;
        for (j = 0; j < r->n; ++j)
            if (r->r[j].wl) { dict_destroy(r->r[j].wl); r->r[j].wl = NULL; }
        free(r->r);
        if (r->raw_tag) free(r->raw_tag);
        if (r->corr_tag) free(r->corr_tag);
    }
    
    free(args.bcs);
    free(args.r1->r);
    free(args.r1);
    
    if (args.r2) {
        free(args.r2->r);
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
        b->flag = FQ_FLAG_PASS;
        int exact = 0;
        for (j = 0; j < args.n_bc; ++j) {
            struct bc_reg *r = &args.bcs[j];
            kstring_t str = {0,0,0};
            kstring_t corr = {0,0,0};
            
            int k;
            for (k = 0; k < r->n; ++k) {
                struct bc_reg0 *r0 = &r->r[k];
                char *val = bseq_subset_seq(b, r0->rd, r0->st, r0->ed);
                if (r->corr_tag) {
                    int ex;
                    char *val0 = correct_bc(r0->wl, val, &ex);
                    if (ex) exact++;
                    // if (ex && b->flag == FQ_FLAG_PASS) b->flag = FQ_FLAG_BC_EXACTMATCH;
                    if (val0 == NULL) {
                        b->flag = FQ_FLAG_BC_FAILURE;
                    } else {
                        kputs(val0, &corr);
                        // free(val0);
                    }
                }
                kputs(val, &str);
                free(val);
            }
            // make sure all segments exact matched
            if (exact == r->n)  b->flag = FQ_FLAG_BC_EXACTMATCH;
            
            char *name1 = fname_update_tag(b->n0.s, r->raw_tag, str.s);
            if (corr.l) {
                char *tmp = name1;
                //debug_print("%s", tmp);
                name1 = fname_update_tag(tmp,r->corr_tag, corr.s);
                free(tmp);
            }
            // if corrected, create new tag, otherwise keep empty
            b->n0.l = 0;
            kputs(name1, &b->n0);
            
            free(name1);
            free(str.s);
            if (corr.m) free(corr.s);
        }

        char *r1 = NULL;
        char *q1 = NULL;
        char *r2 = NULL;
        char *q2 = NULL;
        r1 = bseq_subset_seq(b, args.r1->r->rd, args.r1->r->st, args.r1->r->ed);
        q1 = bseq_subset_qual(b, args.r1->r->rd, args.r1->r->st, args.r1->r->ed);

        if (r1 == NULL) error("Empty read one.");
        if (args.r2) {
            r2 = bseq_subset_seq(b, args.r2->r->rd, args.r2->r->st, args.r2->r->ed);
            q2 = bseq_subset_qual(b, args.r2->r->rd, args.r2->r->st, args.r2->r->ed);
        }
        
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

void fastq_parse_order()
{
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
}

void fastq_parse_unorder()
{
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
}

extern int fastq_parse2_usage();
int fastq_parse2(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return fastq_parse2_usage();

    if (args.order)
        fastq_parse_order();
    else
        fastq_parse_unorder();
    
    write_report();
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
