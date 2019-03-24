#include "utils.h"
#include "thread_pool.h"
#include "htslib/sam.h"
#include "htslib/khash.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "number.h"

KHASH_SET_INIT_STR(name)

static int usage()
{
    fprintf(stderr, "bam_rmdup input.srt.bam -o out.bam\n");
    fprintf(stderr, "   -tag        Sample tag, cell barcode tag, and/or UMI tag. RG,CB,UR\n");
    fprintf(stderr, "   -t          Threads.\n");
    fprintf(stderr, "   -o          Output bam.\n");
    fprintf(stderr, "   -r          Records per thread chunk.[10000000]\n");
    fprintf(stderr, "   -kd         Keep duplicates, make flag instead of remove them.\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *tag_string;
    int n_thread;
    int keep_dup;
    int bufsize;
    int n_tag;
    char **tags;
    
    bam1_t *last_bam;

    htsFile *fp;
    BGZF *out;
    bam_hdr_t *hdr;
    
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .tag_string = NULL,
    .n_thread = 5,
    .keep_dup = 0,
    .bufsize = 1000000, // 10M
    .last_bam = NULL,
    .n_tag = 0,
    .tags = NULL,
    
    .fp = NULL,
    .out = NULL,
    .hdr = NULL,
};
static int parse_args(int argc, char **argv)
{
    int i;
    const char *thread = NULL;
    const char *chunk_size = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag_string;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-r") == 0) var = &chunk_size;
        else if (strcmp(a, "-kd") == 0) {
            args.keep_dup = 1;
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
    if (args.tag_string == NULL) error("No tag specified.");
    
    if (thread) args.n_thread = str2int((char*)thread);
    if (chunk_size) args.bufsize = str2int((char*)chunk_size);

    kstring_t str = {0,0,0};
    kputs(args.tag_string, &str);
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
}
// todo: improve allocated and free SAM pool
struct sam_pool {
    int n, m;
    bam1_t **bam;
    struct args *opts;
};
static void push_sam_pool(bam1_t *b, struct sam_pool *p)
{
    if (p->n == p->m) {
        p->m = p->m == 0 ? args.bufsize : p->m + 100;
        p->bam = realloc(p->bam, p->m*sizeof(void*));
    }
    p->bam[p->n++] = b;
}
// load each read block into buffer, this design may be uneffective for high coverage data, such as wgs
static struct sam_pool *sam_pool_read(htsFile *fp, int bufsize)
{
    struct sam_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    p->opts = &args;
    int i = 0;
    int last_tid = -1;
    int last_end = -1;
    if (args.last_bam) {
        push_sam_pool(args.last_bam, p);
        args.last_bam = NULL; // set point to NULL
    }
    int ret;
    for (;;) {
        bam1_t *b = bam_init1();
        bam1_core_t *c = &b->core;
        ret = sam_read1(fp, args.hdr, b);
        if (ret < 0) break;
        if (last_tid == -1) last_tid = c->tid;
        if (last_tid != c->tid) { // for different chrom, donot put into one pool
            last_end = -1;
            last_tid = -1;
            args.last_bam = b;
            break;
        }
        
        if (c->mpos > c->pos && last_end < c->mpos) last_end = c->mpos;        

        //debug_print("%s\t%d\t%d\t%d", args.hdr->target_name[c->tid],c->pos+1, c->isize, last_end);
        
        if (i >= bufsize) { // start check ends
            if (last_end > 0 && c->tid == last_tid && c->pos > last_end && c->mpos > c->pos) {
                args.last_bam = b;
                break;
            }
        }
        
        push_sam_pool(b, p);
        i++;        
    }
    debug_print("last_end: %d, p->n:%d, bufsize: %d", last_end, p->n, bufsize);
    if (p->n == 0) {
        free(p);
        return NULL;
    }
    return p;
}

static void write_out(struct sam_pool *p)
{
    int i;
    struct args *opts = p->opts;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = p->bam[i];
        if (b->core.flag & BAM_FDUP) {
            if (opts->keep_dup == 1) {
                if (bam_write1(opts->out, b) == -1) error("Failed to write.");
            }
        }
        else {
            if (bam_write1(opts->out, b) == -1) error("Failed to write.");
        }
        bam_destroy1(b);
    }
    free(p->bam);
    free(p);
}
struct sam_stack_buf {
    int n, m;
    bam1_t **p; // point to sam_pool::bam
};

// credit to samtools::bam_rmdup.c
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
    const bam1_core_t *c = &b->core;
    int i;
    for (i = 0; i < args.n_tag; ++i) {
        uint8_t *tag = bam_aux_get(b, args.tags[i]);
        if (!tag) error("No %s tag at alignment. %d:%d", args.tags[i], c->tid, c->pos+1);
        kputs((char*)tag, &str);
    }
    return str.s;
}
static void dump_best(struct sam_stack_buf *buf, khash_t(name) *best_first, struct args *opts)
{
    // pick read id to  best_first
    bam1_t *pp = NULL; // point to best
    int ret;
    khint_t k;
    //  debug_print("stack->n %d, %d", buf->n, buf->p[0]->core.pos +1);
    if (buf->n == 0) error("Empty stack.");
    int i, j = 0;
    int isize = -1;
    char *last_tag = NULL;
    for (;;) {
        // for reads start from same location, but has different isize or barcodes, we divide reads has same isize and barcodes into one group
        // isize = buf->p[j]->core.isize;
        for (i = j; i < buf->n; ++i) {

            if (i==j) j = -1; // reset j in this cycle and set to iter at end of this cycle of next group

            // if record has been check, set point to NULL and skip at next cycle
            if (buf->p[i] == NULL) continue;

            bam1_t *b = buf->p[i];
            bam1_core_t *c = &b->core;
            
            if (c->mpos < c->pos) { // mated reads check already, check best read hash
                k = kh_get(name, best_first, bam_get_qname(b));
                if (k == kh_end(best_first)) {
                    c->flag |= BAM_FDUP;
                }
                buf->p[i] = NULL;
            }
            else {
                char *tag_string = pick_tag_name(b, opts->n_tag, opts->tags);
                if (pp == NULL) {
                    isize = c->isize;
                    pp = b;
                    if (last_tag) free(last_tag);
                    last_tag = tag_string;
                    continue;
                }
                if (c->isize == isize && strcmp(last_tag, tag_string) == 0) {
                    
                    if (sum_qual(pp) > sum_qual(b)) {
                        c->flag |= BAM_FDUP;
                    }
                    else {
                        bam1_core_t *c1 = &pp->core;
                        c1->flag |= BAM_FDUP;
                        pp = b;
                    }
                    buf->p[i] = NULL; // no need check it again
                }
                else {
                    if (j == -1) j = i;
                }
                free(tag_string);
            }
        }
        // debug_print("j: %d, isize: %d",j, isize);
        if (pp) {
            k = kh_put(name, best_first, bam_get_qname(pp), &ret);
            // debug_print("%s", bam_get_qname(pp));
            pp = NULL;
            free(last_tag);
            last_tag = NULL;
        }

        if (j == -1) break; // no more groups
    }

    buf->n = 0; // clear stack
}
static void push_stack(struct sam_stack_buf *buf, bam1_t *b)
{
    if (buf->n == buf->m) {
        buf->m = buf->m == 0 ? 10 : buf->m*2;
        buf->p = realloc(buf->p, buf->m*sizeof(void*));
    }
    buf->p[buf->n++] = b;
}
static void *run_it(void *_p, int idx)
{
    struct sam_pool *p = (struct sam_pool*)_p;
    struct sam_stack_buf buf;
    memset(&buf, 0, sizeof(struct sam_stack_buf));
    struct args *opts = p->opts;
    // debug_print("p->n: %d", p->n);
    khash_t(name) *best_first = kh_init(name);
    
    int i;
    //bam1_t *lbam = NULL;
    int last_tid = -2;
    int last_pos = -1;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = p->bam[i];
        if (b == NULL)error("Try to access empty point.");
        bam1_core_t *c = &b->core;
        if (last_pos == -1) last_pos = c->pos;
        if (last_tid == -1) last_tid = c->tid;
        if (c->tid != last_tid || last_pos != c->pos) {
            if (buf.n > 0) dump_best(&buf, best_first, opts);
            last_tid = c->tid;
            last_pos = c->pos;
        }

        // for unmapped reads, different chromosomes, singletons, append to output
        if (c->tid == -1) break; // append to output, for sorted BAM unmapped reads only happened at end of file
        if ( !(c->flag&BAM_FPAIRED) || (c->flag&(BAM_FUNMAP|BAM_FMUNMAP)) || (c->mtid >= 0 && c->tid != c->mtid)) continue;
        push_stack(&buf, b);
    }
    if (buf.n > 0) dump_best(&buf, best_first, opts);
    kh_destroy(name, best_first);
    return p;
}

int bam_rmdup(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();

    if (args.n_thread == 1) {
        for (;;) {
            struct sam_pool *b = sam_pool_read(args.fp, args.bufsize);
            if (b == NULL) break;
            b = run_it(b, 0);
            write_out(b);
        }
    }
    else {
        struct thread_pool *p = thread_pool_init(args.n_thread);
        struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
        struct thread_pool_result *r;

        for (;;) {
            struct sam_pool *b = sam_pool_read(args.fp, args.bufsize);
            if (b == NULL) break;
            int block;
            do {
                block = thread_pool_dispatch2(p, q, run_it, b, 0);
                if ((r = thread_pool_next_result(q))) {
                    struct sam_pool *d = (struct sam_pool*)r->data;
                    write_out(d);
                }
            }
            while (block == -1);
        }
        thread_pool_process_flush(q);

        while ((r = thread_pool_next_result(q))) {
            struct sam_pool *d = (struct sam_pool*)r->data;
            write_out(d);
        }

        thread_pool_process_destroy(q);
        thread_pool_destroy(p);
    }
    memory_release();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    
    return 0;
}
