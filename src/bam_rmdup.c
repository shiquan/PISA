#include "utils.h"
#include "dict.h"
#include "htslib/thread_pool.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "number.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *report_fname;
    
    int file_thread;
    int keep_dup;

    int n_tag;
    char **tags;
    htsFile *fp;
    BGZF *out;
    FILE *fp_report;
    bam_hdr_t *hdr;
    int as_SE;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .report_fname = NULL,
    .file_thread  = 1,
    .keep_dup     = 0,
    .n_tag        = 0,
    .tags         = NULL,    
    .fp           = NULL,
    .out          = NULL,
    .fp_report    = NULL,
    .hdr          = NULL,
    .as_SE        = 0,
};
static int parse_args(int argc, char **argv)
{
    int i;
    
    const char *tag_str  = NULL;
    const char *file_thread = NULL;
    
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-tag") == 0) var = &tag_str;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-k") == 0) {
            args.keep_dup = 1;
            continue;
        }
        else if (strcmp(a, "-S") == 0) {
            args.as_SE = 1;
            continue;
        }
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument: %s", a);
    }

    if (args.input_fname == NULL) error("No input BAM.");
    if (args.output_fname == NULL) error("No output BAM specified.");
    if (tag_str == NULL) error("No tag specified.");

    if (file_thread) args.file_thread = str2int((char*)file_thread);
    
    kstring_t str = {0,0,0};
    kputs(tag_str, &str);
    int *s = ksplit(&str, ',', &args.n_tag);
    if (args.n_tag == 0) error("No tags.");
    int j;
    args.tags = malloc(args.n_tag*sizeof(char*));
    for (j = 0; j < args.n_tag; ++j)
        args.tags[j] = strdup(str.s+s[j]);
    free(s);
    free(str.s);
    
    args.fp = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.fp, "%s : %s.", args.input_fname, strerror(errno));
    
    htsFormat type = *hts_get_format(args.fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support SAM/BAM/CRAM format.");
    
    args.hdr = sam_hdr_read(args.fp);
    CHECK_EMPTY(args.hdr, "Failed to open header.");

    if (args.file_thread > 1)
        hts_set_threads(args.fp, args.file_thread);

    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        CHECK_EMPTY(args.fp_report, "%s : %s.", args.report_fname, strerror(errno));
    }
    args.out = bgzf_open(args.output_fname, "w");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    if (bam_hdr_write(args.out, args.hdr) == -1) error("Failed to write SAM header.");
    
    return 0;
}
static void memory_release()
{
    hts_close(args.fp);
    bgzf_close(args.out);
    bam_hdr_destroy(args.hdr);
    int i;
    for (i = 0; i < args.n_tag; ++i) free(args.tags[i]);
    free(args.tags);
    if (args.fp_report) fclose(args.fp_report);
}

static int all_reads = 0;
static int duplicate = 0;

static inline int sum_qual(const bam1_t *b)
{
    int i, q;
    uint8_t *qual = bam_get_qual(b);
    for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
    return q;
}
static inline char *pick_tag_name(const bam1_t *b, int n_tag, char **tags)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < args.n_tag; ++i) {
        uint8_t *tag = bam_aux_get(b, args.tags[i]);        
        if (!tag){
            //error("No %s tag at alignment. %d:%lld", args.tags[i], c->tid, c->pos+1);
            if (str.l) free(str.s);
            return NULL;
        }
        kputs((char*)tag, &str);
    }
    return str.s;
}

static struct {
    int n, m;
    bam1_t **b;
    struct dict *best_names;
} buf = {
    .n = 0,
    .m = 0,
    .b = NULL,
    .best_names = NULL,
};

static void clean_buffer()
{
    int i;
    for (i = 0; i < buf.n; ++i) bam_destroy1(buf.b[i]);
    buf.n = 0;
}
static void clean_buffer1()
{
    clean_buffer();
    if (buf.best_names) dict_destroy(buf.best_names);
    buf.best_names = NULL;
}
static void destroy_buffer()
{
    clean_buffer1();
    if (buf.m) free(buf.b);
}
static void push_buffer(bam1_t *b)
{
    if (buf.n == buf.m) {
        buf.m = buf.n == 0 ? 12 : buf.m * 2;
        buf.b = realloc(buf.b, buf.m *sizeof(void*));
    }
    buf.b[buf.n++] = bam_dup1(b);
}

struct read_qual {
    char *name; // point to read struct
    int qual;
};

struct rq_groups {
    struct rq_groups *next;
    int isize;
    int n, m;
    struct read_qual *q;
};

static void dump_best()
{
    if (buf.n == 0) return;
    
    if (buf.n == 1) {
        if (bam_write1(args.out, buf.b[0]) == -1) error("Failed to write.");
        clean_buffer();

        bam1_core_t *c = &buf.b[0]->core;
        if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY)
            return;
        
        all_reads++;
        return;
    }
    //debug_print("%d", buf.n);
    struct dict *reads_group = dict_init();
    dict_set_value(reads_group);
    
    if (buf.best_names == NULL) buf.best_names = dict_init();
    
    int i;

    // groups reads from the same fragment
    for (i = 0; i < buf.n; ++i) {
        bam1_t *b = buf.b[i];
        bam1_core_t *c = &b->core;
        
        if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY) continue;
        
        int isize = c->isize;
        if (isize != 0 && args.as_SE == 1) isize = 0;
        if (isize == 0) { // update isize to read length, for SE mode
            int endpos = bam_endpos(b);
            isize = endpos - c->pos;
        }
        char *bc = pick_tag_name(b, args.n_tag, args.tags);

        if (bc == NULL) continue;
        int idx = dict_query(reads_group, bc);

        if (idx == -1) idx = dict_push(reads_group, bc);
        struct rq_groups *r = dict_query_value(reads_group, idx);
        struct rq_groups *last_r = NULL;
        for (;;) {
            if (r == NULL) {
                r = malloc(sizeof(*r));
                r->n = r->m = 0;
                r->next = NULL;
                r->q = NULL;
                r->isize = isize; // init isize
                if (last_r == NULL) // header
                    dict_assign_value(reads_group, idx, r);
                    
                else 
                    last_r->next = r;
                last_r = r;
            }
            
            if (isize != r->isize) r = r->next;
            else break; 
        }
        if (r->n == r->m) {
            r->m = r->n == 0 ? 2 : r->m *2;
            r->q = realloc(r->q, r->m*sizeof(struct read_qual));
        }
        r->q[r->n].name = bam_get_qname(b);
        r->q[r->n].qual = sum_qual(b);
        
        r->n++;
        free(bc);
    }

    // select read name with best quality
    for (i = 0; i < dict_size(reads_group); ++i) {
        struct rq_groups *r = dict_query_value(reads_group, i);
        assert(r);
        while (r) {

            if (r->isize < 0) {
                r = r->next;
                continue;
            }
            int j;
            int best_read = 0;
            int qual = -1;
            for (j = 0; j < r->n; ++j) {
                if (qual < r->q[j].qual) {
                    qual = r->q[j].qual;
                    best_read = j;
                }
            }
            dict_push(buf.best_names, r->q[best_read].name); // keep best read name
            //debug_print("Best name, %s", r->q[best_read].name);
            r = r->next;
        }
    }

    // export reads
    for (i = 0; i < buf.n; ++i) {
        bam1_t *b = buf.b[i];
        bam1_core_t *c = &b->core;
        if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY) {
            if (bam_write1(args.out, b) == -1) error("Failed to write.");
            continue;
        }

        all_reads++;
        
        int idx = dict_query(buf.best_names, bam_get_qname(b));

        if (idx < 0) {
            c->flag |= BAM_FDUP;
            duplicate++;
        }                
        if (args.keep_dup == 0 && idx < 0) continue;
        if (bam_write1(args.out, b) == -1) error("Failed to write.");
    }

    // free
    for (i = 0; i < dict_size(reads_group); ++i) {
        struct rq_groups *r = dict_query_value(reads_group, i);
        while (r) {
            free(r->q);
            void *r1 = r;
            r = r->next;
            free(r1);
        }
    }
    dict_destroy(reads_group);
    clean_buffer();
}
static void print_unmapped(bam1_t *b)
{
    dump_best();
    if (bam_write1(args.out, b) == -1)
        error("Failed to write.");
}
static void summary_report()
{
    if (args.fp_report) {
        fprintf(args.fp_report, "All reads,%d\n", all_reads);
        fprintf(args.fp_report, "Duplicate reads,%d\n", duplicate);
        fprintf(args.fp_report, "Duplicate ratio,%.4f", (float)duplicate/all_reads);
    }
    LOG_print("All reads,%d", all_reads);
    LOG_print("Duplicate reads,%d", duplicate);
    LOG_print("Duplicate ratio,%.4f", (float)duplicate/all_reads);
}
extern int rmdup_usage();

int bam_rmdup(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return rmdup_usage();

    bam1_t *b = bam_init1();
    bam1_core_t *c = &b->core;
    int ret;
    int last_tid = -2;
    int last_pos = -1;

    for (;;) {
        ret = sam_read1(args.fp, args.hdr, b);
        if (ret < 0) break; // end of file

        // assume inputs are sorted
        if (c->tid == -1) {
            print_unmapped(b);
            continue;
        }
        if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY) {
            push_buffer(b);
            continue;
        }


        if (last_tid != c->tid) {
            LOG_print("Deduplicating %s", args.hdr->target_name[c->tid]);
            last_tid = c->tid;
            last_pos = -1;
            dump_best();
            clean_buffer1();
        }

        if (last_pos == -1) {
            last_pos = c->pos;
        }
        else if (last_pos == c->pos) {

        }
        else if (last_pos < c->pos) {
            dump_best();
        }
        else {
            error("Unsorted bam?");
        }
        last_pos = c->pos;
        push_buffer(b);
    }
    dump_best();
    destroy_buffer();
    
    summary_report();
    
    bam_destroy1(b);
    memory_release();
    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    
    return 0;
}
