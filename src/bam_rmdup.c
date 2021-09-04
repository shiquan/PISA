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
    int qual_thres;
    int n_tag;
    char **tags;
    htsFile *fp;
    BGZF *out;
    FILE *fp_report;
    bam_hdr_t *hdr;
    //int as_SE;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .report_fname = NULL,
    .file_thread  = 1,
    .keep_dup     = 0,
    .qual_thres   = 0,
    .n_tag        = 0,
    .tags         = NULL,    
    .fp           = NULL,
    .out          = NULL,
    .fp_report    = NULL,
    .hdr          = NULL,
    //.as_SE        = 0,
};
static int parse_args(int argc, char **argv)
{
    int i;
    
    const char *tag_str  = NULL;
    const char *file_thread = NULL;
    const char *qual_thres = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-tags") == 0 || strcmp(a, "-tag") == 0) var = &tag_str;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;
        else if (strcmp(a, "-k") == 0) {
            args.keep_dup = 1;
            continue;
        }
        //else if (strcmp(a, "-S") == 0) {
        //    args.as_SE = 1;
        //    continue;
        //}
        
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
    if (qual_thres) args.qual_thres = str2int((char*)qual_thres);
    
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
    const bam1_core_t *c = &b->core;
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < args.n_tag; ++i) {
        uint8_t *tag = bam_aux_get(b, args.tags[i]);        
        if (!tag){
            error("No %s tag at alignment. %d:%lld", args.tags[i], c->tid, c->pos+1);
            //if (str.l) free(str.s);
            //return NULL;
        }
        kputs((char*)tag, &str);
    }
    return str.s;
}

struct rcd {
    bam1_t *b;
    int idx;
    int qual;
    int dup:1;
    int checked:1;
};

static struct {
    int n, m;
    struct rcd  *r;
    struct dict *bcs;
} buf = {
    .n = 0,
    .m = 0,
    .r = NULL,
    .bcs = NULL,
};

static void clean_buffer()
{
    int i;
    for (i = 0; i < buf.n; ++i) {
        bam_destroy1(buf.r[i].b);
        buf.r[i].idx = -1;
        buf.r[i].qual = -1;
        buf.r[i].dup= 0;
        buf.r[i].checked = 0;
    }
    buf.n = 0;
    if (buf.bcs) dict_destroy(buf.bcs);
    buf.bcs = NULL;
}

static void destroy_buffer()
{
    if (buf.m) free(buf.r);
}
static void push_buffer(bam1_t *b)
{
    if (buf.n == buf.m) {
        buf.m = buf.n == 0 ? 12 : buf.m * 2;
        buf.r = realloc(buf.r, buf.m *sizeof(struct rcd));
    }
    struct rcd *r = &buf.r[buf.n];
    
    r->b = bam_dup1(b);
    bam1_core_t *c = &b->core;
    if (c->qual < args.qual_thres) r->idx = -1;
    else if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY)
        r->idx = -1;
    else {
        r->qual = sum_qual(b);
        r->dup = 0;
        
        char *bc = pick_tag_name(b, args.n_tag, args.tags);        
        if (buf.bcs == NULL) buf.bcs = dict_init();
        r->idx = dict_query(buf.bcs, bc);
        if (r->idx == -1) r->idx = dict_push(buf.bcs, bc);
        free(bc);
    }
    buf.n++;
}
static void dump_best()
{
    if (buf.n == 0) return;

    int i, j;
    for (i = 0; i < buf.n; ++i) {
        struct rcd *r1 = &buf.r[i];
        if (r1->idx == -1) continue;
        if (r1->checked == 1) continue;
        int best = i;
        for (j = i+1; j < buf.n; ++j) {
            struct rcd *r2 = &buf.r[j];
            if (r2->idx == -1 || r2->checked == 1) continue;
            if (r2->idx != r1->idx) continue;
            r1 = &buf.r[best];
            if (r2->qual > r1->qual) {
                best = j; // update the best quality
                r1->dup = 1; // mark last best one dup
            } else {
                r2->dup = 1; // mark current one dup
            }
            r2->checked = 1;
        }
    }

    for (i = 0; i < buf.n; ++i) {
        struct rcd *r = &buf.r[i];
        if (r->b->core.qual < args.qual_thres) continue;
        if (r->idx != -1) all_reads++;
        if (r->dup == 1) {
            if (args.keep_dup == 0) continue;
            r->b->core.flag |= BAM_FDUP;
            if (bam_write1(args.out, r->b) == -1) error("Failed to write.");
        } else {
            
            if (bam_write1(args.out, r->b) == -1) error("Failed to write.");
        }
    }
    
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
        fprintf(args.fp_report, "Duplicate ratio,%.4f\n", (float)duplicate/all_reads);
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
    const bam1_core_t *c = &b->core;
    int ret;
    int last_tid = -2;
    int last_pos = -1;

    for (;;) {
        ret = sam_read1(args.fp, args.hdr, b);
        if (ret < 0) break; // end of file

        if (c->qual < args.qual_thres) continue;
        
        // assume inputs are sorted
        if (c->tid == -1) {
            print_unmapped(b);
            continue;
        }

        if (last_tid != c->tid) {
            LOG_print("Deduplicating %s", args.hdr->target_name[c->tid]);
            last_tid = c->tid;
            last_pos = -1;
            dump_best();
        }

        if (last_pos == -1) {
            last_pos = c->pos;
        }
        else if (last_pos == c->pos) {
            // just push to buffer
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
