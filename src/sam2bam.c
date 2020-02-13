// Convert SAM records to BAM and parse barcode tag from read name to SAM attributions
#include "utils.h"
#include "thread_pool.h"
#include "number.h"
#include "barcode_list.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 16384)

// flag to skip
#define FLG_USABLE 0
#define FLG_MITO 1
#define FLG_FLT  2

// summary structure for final report
struct reads_summary {
    uint64_t n_reads, n_mapped, n_pair_map, n_pair_all, n_pair_good;
    uint64_t n_sgltn, n_read1, n_read2;
    uint64_t n_diffchr, n_pstrand, n_mstrand;
    uint64_t n_qual;
    uint64_t n_mito;
    //uint64_t n_usable;
    uint64_t n_failed_to_parse;
};
static struct reads_summary *reads_summary_create()
{
    struct reads_summary *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(*s));
    return s;
}
static struct reads_summary *merge_reads_summary(struct reads_summary **rs, int n)
{
    struct reads_summary *s = reads_summary_create();

    int i;
    for (i = 0; i < n; ++i) {

        struct reads_summary *s0 = rs[i];

        s->n_reads     += s0->n_reads;
        s->n_mapped    += s0->n_mapped;
        s->n_pair_map  += s0->n_pair_map;
        s->n_pair_all  += s0->n_pair_all;
        s->n_pair_good += s0->n_pair_good;
        s->n_sgltn     += s0->n_sgltn;
        s->n_read1     += s0->n_read1;
        s->n_read2     += s0->n_read2;
        s->n_diffchr   += s0->n_diffchr;
        s->n_pstrand   += s0->n_pstrand;
        s->n_mstrand   += s0->n_mstrand;
        s->n_qual      += s0->n_qual;
        s->n_mito      += s0->n_mito;
        //s->n_usable    += s0->n_usable;
    }
    return s;
}

static struct args {
    // file names
    const char *input_fname;  // alignment input, only SAM format is required
    const char *output_fname; // BAM output, only support BAM format, default is stdout
    //const char *filter_fname; // if filter BAM file is set, filtered reads will output in this file
    const char *report_fname; // summary report for whole file

    const char *mito; // mitochrondria name, the mito ratio will export in the summary report
    const char *mito_fname; // if set, mitochrondria reads passed QC will be exported in this file
    
    int n_thread;
    int buffer_size;  // buffered records in each chunk
    int file_th;
    gzFile fp;        // input file handler
    kstream_t *ks;    // input streaming
    htsFile *fp_out;     // output file handler
    //BGZF *fp_filter;  // filtered reads file handler
    BGZF *fp_mito;    // if not set, mito reads will be treat at filtered reads
    FILE *fp_report;  // report file handler
    
    bam_hdr_t *hdr;   // bam header structure of input

    char *preload_record;

    struct reads_summary **thread_data; // summary data for each thread

    int fixmate_flag; // if set, fix the mate relationship of paired reads first

    int qual_thres; // mapping quality threshold to filter

    int mito_id;
    int PE_flag;
//    int keep_all;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    //.filter_fname = NULL,
    .report_fname = NULL,
    //.btable_fname = NULL,
    .mito = "chrM",
    .mito_fname = NULL,
    
    .n_thread = 1,
   .buffer_size = 1000000, // 1M
    .file_th  = 1,
    .fp = NULL,
    .ks = NULL,
    .fp_out = NULL,
    //.fp_filter = NULL,
    .fp_mito = NULL,
    .fp_report = NULL,
    
    .hdr = NULL,
    .preload_record = NULL,
    
    .thread_data = NULL,

    .fixmate_flag = 0,
    .qual_thres = 10,
    .mito_id = -2,
    .PE_flag = 0,
    //  .keep_all = 0,
};

// Buffer input and output records in a memory pool per thread
struct sam_pool {
    struct args *opts; // point to args
    int n;
    kstring_t **str;
    bam1_t **bam; // bam structure
    int *flag; // export flag
};

static struct sam_pool* sam_pool_init(int buffer_size)
{
    struct sam_pool *p = malloc(sizeof(*p));
    p->n = 0;

    p->str  = malloc(buffer_size*sizeof(void*));
    p->bam  = malloc(buffer_size*sizeof(void*));
    p->flag = malloc(buffer_size*sizeof(int));
    
    memset(p->str,  0, buffer_size*sizeof(void*));
    memset(p->bam,  0, buffer_size*sizeof(void*));
    memset(p->flag, 0, buffer_size*sizeof(int));
    
    return p;
}
static void sam_pool_destroy(struct sam_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        if (p->str[i]) {
            if (p->str[i]->m) free(p->str[i]->s);
            free(p->str[i]);
        }
        if (p->bam[i]) bam_destroy1(p->bam[i]);        
    }
    free(p->str);
    free(p->bam);
    free(p->flag);
    free(p);
}
static struct sam_pool* sam_pool_read(kstream_t *s, int buffer_size)
{
    struct sam_pool *p = sam_pool_init(buffer_size);
    
    kstring_t str = {0,0,0};
    int ret;
    if (args.preload_record) {
        kstring_t *t = malloc(sizeof(*t));
        memset(t, 0, sizeof(kstring_t));
        kputs(args.preload_record, t);
        p->str[0] = t;
        p->bam[0] = bam_init1();
        p->n++;
        free(args.preload_record);
        args.preload_record= NULL;
    }

    for (;;) {
        if (ks_getuntil(s, 2, &str, &ret) < 0) break;
        if (p->n >= buffer_size) { // in case check paired reads name
            // check the read name
            kstring_t *s = p->str[p->n-1];
            int _i;
            for (_i = 0; _i < str.l; ++_i)
                if (str.s[_i] == '|' || isspace(str.s[_i])) break;
            if (strncmp(str.s, s->s, _i) != 0) {
                args.preload_record = strndup(str.s, str.l);
                break;
            }
            p->str = realloc(p->str, sizeof(void*)*(p->n+1));
            p->bam = realloc(p->bam, sizeof(void*)*(p->n+1));
            p->flag = realloc(p->flag, sizeof(int)*(p->n+1));
        }
        

        // skip header
        if (str.s[0] == '@') continue;

        kstring_t *t = malloc(sizeof(*t));
        t->l = str.l;
        t->m = str.m;
        t->s = strdup(str.s);
        p->str[p->n] = t;
        // init for output
        p->bam[p->n] = bam_init1();

        p->n++;
    }

    free(str.s);
    if (p->n == 0) {
        sam_pool_destroy(p);
        return NULL;
    }
    
    return p;
}
static bam_hdr_t *sam_parse_header(kstream_t *s, kstring_t *line)
{
    bam_hdr_t *h = NULL;
    kstring_t str = {0,0,0};
    int ret;
    
    while (ks_getuntil(s, 2, line, &ret) >= 0) {
        if (line->s[0] != '@') break;
        kputsn(line->s, line->l, &str);
        kputc('\n', &str);
    }

    if (ret < -1) goto failed_parse;
    if (str.l == 0) kputsn("", 0, &str);
    h = sam_hdr_parse(str.l, str.s);
    h->l_text = str.l;
    h->text = str.s;
    return h;

  failed_parse:
    bam_hdr_destroy(h);
    free(str.s);
    
    return NULL;
}
static void write_out(struct sam_pool *p)
{
    int i;
    struct args *opts = p->opts;    
    for (i = 0; i < p->n; ++i) {
        if (p->bam[i] == NULL) continue;

        if (p->flag[i] == FLG_FLT) continue; // filter this alignment for low map quality
        
        if (p->flag[i] == FLG_MITO && opts->fp_mito != NULL) {
            if (bam_write1(opts->fp_mito, p->bam[i]) == -1) error("Failed to write.");
        }
        else {
            if (sam_write1(opts->fp_out, opts->hdr, p->bam[i]) == -1) error("Failed to write.");
        }
    }
    sam_pool_destroy(p);
}
static void summary_report(struct args *opts)
{
    struct reads_summary *summary = merge_reads_summary(opts->thread_data, opts->n_thread);
    if (opts->fp_report) {
        fprintf(opts->fp_report,"Raw reads,%"PRIu64"\n", summary->n_reads);
        fprintf(opts->fp_report,"Mapped reads,%"PRIu64" (%.2f%%)\n", summary->n_mapped, (float)summary->n_mapped/summary->n_reads*100);
        fprintf(opts->fp_report,"Mapped reads (paired),%"PRIu64"\n", summary->n_pair_map);
        fprintf(opts->fp_report,"Properly paired reads,%"PRIu64"\n", summary->n_pair_good);
        fprintf(opts->fp_report,"Singleton reads,%"PRIu64"\n", summary->n_sgltn);
        fprintf(opts->fp_report,"Read 1,%"PRIu64"\n", summary->n_read1);
        fprintf(opts->fp_report,"Read 2,%"PRIu64"\n", summary->n_read2);
        fprintf(opts->fp_report,"Paired reads map on diff chr,%"PRIu64"\n", summary->n_diffchr);
        fprintf(opts->fp_report,"Plus strand,%"PRIu64"\n", summary->n_pstrand);
        fprintf(opts->fp_report,"Minus strand,%"PRIu64"\n", summary->n_mstrand);
        fprintf(opts->fp_report,"Mapping quals above %d,%"PRIu64"\n", opts->qual_thres, summary->n_qual);
        if (opts->mito_id != -1)
            fprintf(opts->fp_report,"Mitochondria ratio,%.2f%%\n", (float)summary->n_mito/summary->n_mapped*100);
        // fprintf(opts->fp_report,"Usable reads (ratio),%"PRIu64" (%.2f%%)\n", summary->n_usable, (float)summary->n_usable/summary->n_reads*100);
        if (summary->n_failed_to_parse > 0)
            fprintf(opts->fp_report,"Failed to parse reads,%"PRIu64"\n", summary->n_failed_to_parse);
    }

    free(summary);
}

static int parse_name_str(kstring_t *s)
{
    // CL100053545L1C001R001_2|||BC:Z:TTTCATGA|||CR:Z:TANTGGTAGCCACTAT|||PL:i:20
    // CL100053545L1C001R001_2 .. CR:Z:TANTGGTAGCCACTAT ..
    int n, i;
    for (n = 0; n < s->l && !isspace(s->s[n]); ++n);
    for (i = 0; i < n && s->s[i] != '|'; ++i);
    
    // int is_name = 1; // first key is name, next is value

    if (i < n-5) { // detected
        char *p = s->s+i;
        char *r = 0;
        char *e = s->s+n;
        kstring_t t = {0,0,0};
        kputsn(s->s, i, &t);
        kputsn(s->s+n, s->l-n, &t);
        *e = '\0';

        for ( ; p != e; ) {
            if (*p == '|' && *(p+1) == '|' && *(p+2) == '|') {
                *p = '\0';                
                if (r != 0) {
                    kputc('\t', &t);
                    kputs(r, &t);
                }
                p += 3;
                r = p;                
            }
            p++;
        }
        if (r) {
            kputc('\t', &t);
            kputs(r, &t);
        }
        s->l = t.l;
        s->m = t.m;
        free(s->s);
        s->s = t.s;
    }

    return 0;
}
static void sam_stat_reads(bam1_t *b, struct reads_summary *s, int *flag, struct args *opts)
{
    bam1_core_t *c = &b->core;
    s->n_reads++;

    if (c->qual >= opts->qual_thres) s->n_qual++;
    else {
        *flag = FLG_FLT; // filter low mapping quality
        return;
    }
            
    if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY || c->flag & BAM_FUNMAP) return; // unmapped


    // mito ratio, mito come before qc filter, so only high quality mito reads will be exported in mito.bam
    if (c->tid == opts->mito_id) {
        s->n_mito++;
        *flag = FLG_MITO; // filter mito
    }
        
    s->n_mapped++;
    if (c->flag & BAM_FREVERSE) s->n_mstrand++;
    else s->n_pstrand++;            

    if (c->flag & BAM_FPAIRED) {
        s->n_pair_all ++;

        if ((c->flag & BAM_FPROPER_PAIR) && !(c->flag & BAM_FUNMAP)) s->n_pair_good++;
            
        if (c->flag & BAM_FREAD1) s->n_read1++;
        else if (c->flag & BAM_FREAD2) s->n_read2++;

        if ((c->flag & BAM_FMUNMAP) && !(c->flag & BAM_FUNMAP)) {
            s->n_sgltn++;
        }
        
        if (!(c->flag & BAM_FUNMAP) && !(c->flag & BAM_FMUNMAP)) {
            s->n_pair_map++;
            if (c->mtid != c->tid) {
                s->n_diffchr++;
                // *flag = FLG_FLT; // no matter -p set or not, skip reads mapped to diff chroms
            }
        }    
    }
    //if (*flag == FLG_USABLE) s->n_usable++;
}
static void *sam_name_parse(void *_p, int idx)
{
    struct sam_pool *p = (struct sam_pool*)_p;
    struct args *opts = p->opts;
    struct reads_summary *summary = opts->thread_data[idx];    
    bam_hdr_t *h = opts->hdr;
    if (opts->PE_flag == 0) {
        int i;
        for (i = 0; i < p->n; ++i) {
            parse_name_str(p->str[i]);
            if (sam_parse1(p->str[i], h, p->bam[i])) {
                warnings ("Failed to parse SAM., %s", bam_get_qname(p->bam[i]));
                summary->n_failed_to_parse++;
                p->bam[i] = NULL;
                continue;
            }
            sam_stat_reads(p->bam[i], summary, &p->flag[i], opts);
        }
    }
    else { // PE
        
        bam1_t *b1 = NULL, *b2 = NULL;
        //int *f1 = NULL, *f2 = NULL;
        int i;
        for (i = 0; i < p->n; ++i) {
            parse_name_str(p->str[i]);
            if (b1 == NULL) {
                b1 = p->bam[i];
                if (sam_parse1(p->str[i], h, b1)) {
                    warnings ("Failed to parse SAM., %s", bam_get_qname(p->bam[i]));
                    summary->n_failed_to_parse++;
                    p->bam[i] = NULL;
                    continue;
                }
                if (b1->core.flag&BAM_FSECONDARY || b1->core.flag&BAM_FSUPPLEMENTARY) {
                    p->flag[i] = FLG_FLT;
                    b1 = NULL;
                    continue;
                }
                sam_stat_reads(b1, summary, &p->flag[i], opts);
                //f1 = &p->flag[i];
            }
            else {
                b2 = p->bam[i];
                if (sam_parse1(p->str[i], h, b2)) {
                    warnings ("Failed to parse SAM., %s", bam_get_qname(p->bam[i]));
                    summary->n_failed_to_parse++;
                    p->bam[i] = NULL;
                    continue;
                    // error("Failed to parse SAM.");
                }
                if (b2->core.flag&BAM_FSECONDARY || b2->core.flag&BAM_FSUPPLEMENTARY) {
                    p->flag[i] = FLG_FLT;
                    b2 = NULL;
                    continue;
                }
                if (strcmp(bam_get_qname(b1), bam_get_qname(b2)) != 0) error("Inconsist paried read name. %s vs %s", bam_get_qname(b1), bam_get_qname(b2));
                sam_stat_reads(b2, summary, &p->flag[i], opts);
                //f2 = &p->flag[i];
                /*
                if (*f1 == FLG_FLT || *f2 == FLG_FLT || b1->core.tid != b2->core.tid) {
                    *f1 = FLG_FLT;
                    *f2 = FLG_FLT;
                }
                */
                b1 = NULL;
                b2 = NULL;
            }     
        }
    }
    
    return p;
}
static int sam_name_parse_light()
{
    for (;;) {
        struct sam_pool *p = sam_pool_read(args.ks, args.buffer_size);
        if (p == NULL) break;
        p->opts = &args;
        p = sam_name_parse(p, 0);
        write_out(p);
    }
    return 0;
}

extern int sam2bam_usage();

static int parse_args(int argc, char **argv)
{
    int i;
    const char *buffer_size = NULL;
    const char *thread = NULL;
    const char *qual_thres = NULL;
    const char *file_th = NULL;
    
    for (i = 1; i < argc;) {

        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return sam2bam_usage();
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-r") == 0) var = &buffer_size;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;
        else if (strcmp(a, "-mito") == 0) var = &args.mito;
        else if (strcmp(a, "-maln") == 0) var = &args.mito_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-filter") == 0) continue; // var = &args.filter_fname;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-k") == 0) { // -k has been removed, 2020/02/13
            // args.keep_all = 1;
            continue; 
        }
        else if (strcmp(a, "-p") == 0) {
            args.PE_flag = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1]) error("Unknown parameter: %s", a);

        if (args.input_fname == 0) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument: %s.", a);
    }

    // init input
    if (args.input_fname == NULL && !isatty(fileno(stdin))) args.input_fname = "-";
    if (args.input_fname == NULL) error("No input SAM file is set!");
    if (args.output_fname == NULL) error("No output BAM file specified.");
    args.fp = strcmp(args.input_fname, "-") ? gzopen(args.input_fname, "r") : gzdopen(fileno(stdin), "r");
    if (args.fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    args.ks = ks_init(args.fp);

    // init output
    
    args.fp_out = hts_open(args.output_fname, "bw");
    if (args.fp_out == NULL) error("%s : %s.", args.output_fname, strerror(errno));

    if (file_th) {
        args.file_th = str2int((char*)file_th);
        if (args.file_th <1) args.file_th = 1;
        hts_set_threads(args.fp_out, args.file_th);
    }
    
    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        if (args.fp_report == NULL) error("%s : %s.", args.report_fname, strerror(errno));
    }
    /*
    if (args.filter_fname) {
        if (args.keep_all == 1) error("-k is conflict with -filter");
        args.fp_filter = bgzf_open(args.filter_fname, "w");
        if (args.fp_filter == NULL) error("%s : %s.", args.filter_fname, strerror(errno));
    }
    */
    if (args.mito_fname) {
        // if (args.keep_all == 1) error("-k is conflict with -maln");
        args.fp_mito = bgzf_open(args.mito_fname, "w");
        if (args.fp_mito == NULL) error("%s : %s.", args.mito_fname, strerror(errno));
    }
    // init parameters

    // the bottleneck happens at compress BAM, not processing sam records, so we delete the multi-thread mode, 2019/11/14
    /*
    if (thread) {
        args.n_thread = str2int((char*)thread);
        assert(args.n_thread > 0);
    }
    */
    if (buffer_size) {
        args.buffer_size = str2int((char*)buffer_size);
        assert(args.buffer_size>0);
    }
    if (qual_thres) {
        args.qual_thres = str2int((char*)qual_thres);
        assert(args.qual_thres >= 0);
    }

    // init thread data
    args.thread_data = malloc(args.n_thread*sizeof(struct reads_summary*));
    for (i = 0; i < args.n_thread; ++i) args.thread_data[i] = reads_summary_create();
        
    // init bam header and first bam record
    kstring_t str = {0,0,0}; // cache first record
    args.hdr = sam_parse_header(args.ks, &str);
    if (args.hdr == NULL) error("Failed to parse header. %s", args.input_fname);
    if (sam_hdr_write(args.fp_out, args.hdr)) error("Failed to write header.");
    // if (args.fp_filter && bam_hdr_write(args.fp_filter, args.hdr)) error("Failed to write header.");
    if (args.fp_mito && bam_hdr_write(args.fp_mito, args.hdr)) error("Failed to write header.");

    // init mitochondria id
    args.mito_id = bam_name2id(args.hdr, args.mito);
    if (args.mito_id == -1) {
        warnings("Failed to find mitochondria name at the header, will NOT stat the mito ratio!");
        args.mito_id = -2; // reset to -2, in case inconsistance with unmap
    }

    // check if there is a BAM record
    if (str.s[0] != '@') {        
        args.preload_record = strndup(str.s, str.l);
        /*
          parse_name_str(&str);        
        bam1_t *b = bam_init1();
        if (sam_parse1(&str, args.hdr, b)) error("Failed to parse SAM.");

        int flag = 0;        
        sam_stat_reads(b, args.thread_data[0], &flag, &args);
        
        if (flag == FLG_USABLE) {
            if (bam_write1(args.fp_out, b) == -1) error("Failed to write BAM.");
        }
        else if (flag == FLG_MITO && args.fp_mito != NULL) {
            if (bam_write1(args.fp_mito, b) == -1) error("Failed to write BAM.");
        }
        else if (args.fp_filter != NULL) {
            if (bam_write1(args.fp_filter, b) == -1) error("Failed to write BAM.");
        }
        else {

        }
        bam_destroy1(b);
        */
    }
    free(str.s);
    return 0;
}

static void memory_release()
{
    hts_close(args.fp_out);
    ks_destroy(args.ks);
    gzclose(args.fp);
    bam_hdr_destroy(args.hdr);

    //if (args.fp_filter) bgzf_close(args.fp_filter);
    if (args.fp_mito) bgzf_close(args.fp_mito);
    if (args.fp_report) fclose(args.fp_report);
    int i;
    for (i = 0; i < args.n_thread; ++i) free(args.thread_data[i]);
    free(args.thread_data);
}

int sam2bam(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return 1;

    if (args.n_thread == 1)
        sam_name_parse_light();
    
    else {

        int nt = args.n_thread;

        struct thread_pool *p = thread_pool_init(nt);
        struct thread_pool_process *q = thread_pool_process_init(p, nt*2, 0);
        struct thread_pool_result *r;

        for (;;) {

            struct sam_pool *b = sam_pool_read(args.ks, args.buffer_size);
            if (b == NULL) break;
            b->opts = &args;
            
            int block;

            do {

                block = thread_pool_dispatch2(p, q, sam_name_parse, b, 0);

                if ((r = thread_pool_next_result(q))) {
                    struct sam_pool *d = (struct sam_pool*)r->data;
                    write_out(d);
                }
                //thread_pool_delete_result(r, 1);
            }
            while (block == -1);
        }
        
        thread_pool_process_flush(q);

        while ((r = thread_pool_next_result(q))) {
            struct sam_pool *d = (struct sam_pool *)r->data;
            write_out(d);
        }

        thread_pool_process_destroy(q);
        thread_pool_destroy(p);
        
    }

    summary_report(&args);
    
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
