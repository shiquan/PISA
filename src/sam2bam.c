// Convert SAM records to BAM and parse barcode tag from read name to SAM attributions
#include "utils.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "thread_pool_internal.h"
#include <zlib.h>
#include "gtf.h"
#include "read_anno.h"

static char *corr_tag = "MM";

KSTREAM_INIT(gzFile, gzread, 16384)

// flag to skip
#define FLG_USABLE 0
#define FLG_MITO 1
#define FLG_FLT  2

// summary structure for final report
struct summary {
    uint64_t n_reads, n_mapped, n_pair_map, n_pair_all, n_pair_good;
    uint64_t n_sgltn, n_read1, n_read2;
    uint64_t n_diffchr, n_pstrand, n_mstrand;
    uint64_t n_mito;
    uint64_t n_failed_to_parse;
    uint64_t n_adj;
    uint64_t n_corr; // mapq corrected
};
static struct summary *summary_create()
{
    struct summary *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(*s));
    return s;
}

static struct args {
    // file names
    const char *input_fname;  // alignment input, only SAM format is required
    const char *output_fname; // BAM output, only support BAM format, default is stdout
    const char *report_fname; // summary report for whole file

    const char *mito; // mitochrondria name, the mito ratio will export in the summary report
    const char *mito_fname; // if set, mitochrondria reads passed QC will be exported in this file

    const char *gtf_fname; // gtf is required if -adjust-mapq set

    int qual_corr;
    int enable_corr;
    struct gtf_spec *G;
    
    int n_thread;
    int buffer_size;  // buffered records in each chunk
    int file_th;
    gzFile fp;        // input file handler
    kstream_t *ks;    // input streaming
    htsFile *fp_out;     // output file handler

    BGZF *fp_mito;    // if not set, mito reads will be treat at filtered reads
    FILE *fp_report;  // report file handler
    
    bam_hdr_t *hdr;   // bam header structure of input

    char *preload_record;

    struct summary *summary;

    int mito_id;
    int qual_thres;
} args = {
    .input_fname       = NULL,
    .output_fname      = NULL,
    .report_fname      = NULL,
    .gtf_fname         = NULL,
    .mito              = "chrM",
    .mito_fname        = NULL,
    .qual_corr         = 255,
    .enable_corr       = 0,
    .G                 = NULL,
    .n_thread          = 1,
    .buffer_size       = 1000000, // 1M
    .file_th           = 1,
    .fp                = NULL,
    .ks                = NULL,
    .fp_out            = NULL,
    .fp_mito           = NULL,
    .fp_report         = NULL,
    .hdr               = NULL,
    .preload_record    = NULL,
    .summary           = NULL,
    .mito_id           = -2,
    .qual_thres        = 0,
};

pthread_mutex_t global_data_mutex = PTHREAD_MUTEX_INITIALIZER;

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

        /* do NOT filter any records, edited 2020/04/04
        if (p->flag[i] == FLG_FLT) continue; // filter this alignment for low map quality
        */
        if (p->flag[i] == FLG_MITO && opts->fp_mito != NULL) {
            if (bam_write1(opts->fp_mito, p->bam[i]) == -1) error("Failed to write.");
            continue;
        }
        if (sam_write1(opts->fp_out, opts->hdr, p->bam[i]) == -1) error("Failed to write.");
    }
    sam_pool_destroy(p);
}
static void summary_report(struct args *opts)
{
    struct summary *summary = opts->summary;
    fprintf(opts->fp_report,"Raw reads,%"PRIu64"\n", summary->n_reads);
    fprintf(opts->fp_report,"Mapped reads,%"PRIu64" (%.2f%%)\n", summary->n_mapped, (float)summary->n_mapped/summary->n_reads*100);
    if (summary->n_pair_map > 0) {
        fprintf(opts->fp_report,"Mapped reads (paired),%"PRIu64"\n", summary->n_pair_map);
        fprintf(opts->fp_report,"Properly paired reads,%"PRIu64"\n", summary->n_pair_good);
        fprintf(opts->fp_report,"Singleton reads,%"PRIu64"\n", summary->n_sgltn);
        fprintf(opts->fp_report,"Read 1,%"PRIu64"\n", summary->n_read1);
        fprintf(opts->fp_report,"Read 2,%"PRIu64"\n", summary->n_read2);
        fprintf(opts->fp_report,"Paired reads map on diff chr,%"PRIu64"\n", summary->n_diffchr);
    }
    fprintf(opts->fp_report,"Plus strand,%"PRIu64"\n", summary->n_pstrand);
    fprintf(opts->fp_report,"Minus strand,%"PRIu64"\n", summary->n_mstrand);
    if (opts->mito_id != -1)
        fprintf(opts->fp_report,"Mitochondria ratio,%.2f%%\n", (float)summary->n_mito/summary->n_mapped*100);
    if (summary->n_failed_to_parse > 0)
        fprintf(opts->fp_report,"Failed to parse reads,%"PRIu64"\n", summary->n_failed_to_parse);
    if (opts->enable_corr)
        fprintf(opts->fp_report, "Mapping quality corrected reads,%"PRIu64"\n", summary->n_corr);
}

static int parse_name_str(kstring_t *s)
{
    // CL100053545L1C001R001_2|||BC:Z:TTTCATGA|||CR:Z:TANTGGTAGCCACTAT|||PL:i:20
    // CL100053545L1C001R001_2 .. CR:Z:TANTGGTAGCCACTAT ..
    int n, i;
    for (n = 0; n < s->l && !isspace(s->s[n]); ++n);
    for (i = 0; i < n && s->s[i] != '|'; ++i);
    
    // int is_name = 1; // first key is name, next is value

    if (i> 0 && i < n-5) { // detected
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
static void sam_stat_reads(bam1_t *b, struct summary *s, int *flag, struct args *opts)
{
    bam1_core_t *c = &b->core;
    if (c->flag & BAM_FSECONDARY) return; // skip secondary alignment
    s->n_reads++;

    if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSUPPLEMENTARY || c->flag & BAM_FUNMAP) return; // unmapped


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
            }
        }
    }    
}
extern struct gtf_anno_type *bam_gtf_anno_core(bam1_t *b, struct gtf_spec const *G, bam_hdr_t *h);
extern void gtf_anno_destroy(struct gtf_anno_type *ann);
extern int sam_realloc_bam_data(bam1_t *b, size_t desired);
// return 0 on not correct, 1 on corrected
static void shrink_bam(bam1_t *bam)
{
    bam1_core_t *c = &bam->core;
    if (c->flag & BAM_FSECONDARY) {
        c->qual = 0;
        if (c->l_qseq > 0) {
            uint8_t *s = bam->data + (c->n_cigar<<2) + c->l_qname;
            uint8_t *e = bam->data + bam->l_data;
            int l_data = ((c->l_qseq+1)>>1)+c->l_qseq; 
            uint8_t *r = s + l_data;
            memmove(s, r, e-r);
            bam->l_data = bam->l_data - l_data;
            c->l_qseq = 0;
        }
    }
}
int bam_map_qual_corr(bam1_t **b, int n, struct gtf_spec const *G, int qual)
{
    int i;
    int best_hits = 0;
    int best_bam  = -1;
    uint8_t *data = NULL;
    int l_data = 0;
    int l_qseq = 0;
    for (i = 0; i < n; ++i) {
        bam1_t *bam = b[i];
        bam1_core_t *c = &bam->core;
        if (l_data == 0 && c->l_qseq > 0) {
            l_data = ((c->l_qseq+1)>>1)+c->l_qseq; 
            data = malloc(l_data);
            memcpy(data, bam->data + (c->n_cigar<<2) + c->l_qname, l_data);
            l_qseq = c->l_qseq;
        }
        struct gtf_anno_type *ann = bam_gtf_anno_core(bam, G, args.hdr);
        if (ann == NULL) continue;
        // read mapped in exon will be selected
        if (ann->type != type_exon &&
            ann->type != type_splice &&
            ann->type != type_exon_intron) {
            gtf_anno_destroy(ann);

            continue;
        }
            
        if (c->flag & BAM_FSECONDARY) best_bam = i;
        best_hits++;
        gtf_anno_destroy(ann);
    }
    // only one secondary alignment hit exonic region
    if (best_hits > 1) {
        if (data) free(data);
        for (i = 0; i < n; ++i) shrink_bam(b[i]);
        return 0;
    }
    if (best_bam == -1) {
        if (data) free(data);
        for (i = 0; i < n; ++i) shrink_bam(b[i]);
        return 0;
    }
    // update the mapping quality and flag
    for (i = 0; i < n; ++i) {
        bam1_t *bam = b[i];
        bam1_core_t *c = &bam->core;
        if (i == best_bam) {
            int flag = BAM_FSECONDARY;            
            c->flag &= ~flag;
            c->qual = qual;
            if (c->l_qseq == 0 && l_data > 0) {
                c->l_qseq = l_qseq;
                int m_data = bam->l_data + l_data;
                if (m_data > bam->m_data) {
                    sam_realloc_bam_data(bam, m_data);
                }
                uint8_t *s = bam->data + (c->n_cigar<<2) + c->l_qname;
                uint8_t *e = bam->data + bam->l_data;
                memmove(s+l_data, s, e-s);
                memcpy(s, data, l_data);
                bam->l_data = m_data;
            }
            bam_aux_update_int(bam, corr_tag, 1);
        }
        else {
            c->flag |= BAM_FSECONDARY;
            shrink_bam(bam);
        }
    }
    free(data);
    return 1;
}
// return 0 on same read, 1 on differnt name
int bam_same(bam1_t *a, bam1_t *b)
{
    if (strcmp((char*)a->data, (char*)b->data) != 0) return 1;
    if ((a->core.flag & BAM_FREAD1) && (b->core.flag & BAM_FREAD2)) return 1;
    if ((a->core.flag & BAM_FREAD2) && (b->core.flag & BAM_FREAD1)) return 1;
    return 0;
}
int bam_pool_qual_corr(struct sam_pool *p)
{
    int i;
    int corred = 0;
    for (i = 0; i < p->n; ) {
        bam1_t *bam = p->bam[i];
        if (bam == NULL) {
            i++;continue;
        }
        
        bam1_core_t *c = &bam->core;
        if (c->tid <= -1 || c->tid > args.hdr->n_targets || (c->flag & BAM_FUNMAP)) {
            i++; continue;
        }
        // only correct multi-hits        
        if (!(c->flag & BAM_FSECONDARY)) {
            i++; continue;
        }

        int st,ed;
        for (st = i-1; st >= 0; --st) 
            if (bam_same(bam, p->bam[st]) != 0) break;
       
        st++; // move point back to same record
        
        for (ed = i+1; ed < p->n; ++ed)
            if (bam_same(bam, p->bam[ed]) != 0) break;

        ed--; // move point back

        i = ed+1; // jump to next read
        
        if (ed-st==0) break; // only one record, no need to corr
        
        int n = ed-st+1;
        bam1_t **b = malloc(n*sizeof(void*));
        
        int j;
        for (j = 0; j < ed-st+1; ++j) b[j] = p->bam[st+j];
        corred += bam_map_qual_corr(b, n, args.G, args.qual_corr);
        free(b); // free stack
    }
    return corred;
}
static int sam_safe_check(kstring_t *str)
{
    int i;
    int s = 0;
    for (i = 0; i < str->l; ++i)
        if (isspace(str->s[i])) s++;
    if (s < 11) return 1; // we need at least 11 columns for SAM
    return 0;
}
static void *sam_name_parse(void *_p)
{
    struct sam_pool *p = (struct sam_pool*)_p;
    struct args *opts = p->opts;
    struct summary *s0 = summary_create();
    bam_hdr_t *h = opts->hdr;

    int i;
    for (i = 0; i < p->n; ++i) {
        parse_name_str(p->str[i]);
        if (sam_safe_check(p->str[i])) {
            warnings("Failed to parse %s", p->str[i]->s);
            s0->n_failed_to_parse++;
            p->bam[i] = NULL;
            continue;
        }
        if (sam_parse1(p->str[i], h, p->bam[i])) {
            warnings ("Failed to parse SAM., %s", bam_get_qname(p->bam[i]));
            s0->n_failed_to_parse++;
            p->bam[i] = NULL;
        }
    }
    int n_corr = 0;
    if (args.enable_corr) 
        n_corr = bam_pool_qual_corr(p);
    
    for (i = 0; i < p->n; ++i) 
        sam_stat_reads(p->bam[i], s0, &p->flag[i], opts);

    pthread_mutex_lock(&global_data_mutex);
    struct summary *s = opts->summary;
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
    s->n_mito      += s0->n_mito;
    s->n_adj       += s0->n_adj;
    s->n_corr      += n_corr;
    pthread_mutex_unlock(&global_data_mutex);
    free(s0);
    return p;
}
static int sam_name_parse_light()
{
    for (;;) {
        struct sam_pool *p = sam_pool_read(args.ks, args.buffer_size);
        if (p == NULL) break;
        p->opts = &args;
        p = sam_name_parse(p);
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
    const char *qual_corr = NULL;
    const char *file_th = NULL;
    const char *qual_thres = NULL;
    for (i = 1; i < argc;) {

        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return sam2bam_usage();
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-r") == 0) var = &buffer_size;
        else if (strcmp(a, "-mito") == 0) var = &args.mito;
        else if (strcmp(a, "-maln") == 0) var = &args.mito_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-qual") == 0) var = &qual_corr;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;
        else if (strcmp(a, "-k") == 0) { // -k has been removed, 2020/02/13
            continue; 
        }
        else if (strcmp(a, "-adjust-mapq") == 0) {
            args.enable_corr = 1;
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

    if (qual_thres)
        args.qual_thres = str2int(qual_thres);
    
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

    if (args.enable_corr) {
        if (args.gtf_fname == NULL) error("-gtf is required if mapping quality correction enabled.");
        args.G = gtf_read_lite(args.gtf_fname);
        if (args.G == NULL) error("GTF is empty.");
    }
    
    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        if (args.fp_report == NULL) error("%s : %s.", args.report_fname, strerror(errno));
    }
    else args.fp_report = stdout;
    
    if (args.mito_fname) {
        args.fp_mito = bgzf_open(args.mito_fname, "w");
        if (args.fp_mito == NULL) error("%s : %s.", args.mito_fname, strerror(errno));
    }

    if (thread) {
        args.n_thread = str2int((char*)thread);
        assert(args.n_thread > 0);
    }

    if (buffer_size) {
        args.buffer_size = str2int((char*)buffer_size);
        assert(args.buffer_size>0);
    }
    if (qual_corr) {
        args.qual_corr = str2int((char*)qual_corr);
        assert(args.qual_corr >= 0);
    }

    args.summary = summary_create();
    
    // init bam header and first bam record
    kstring_t str = {0,0,0}; // cache first record
    args.hdr = sam_parse_header(args.ks, &str);
    if (args.hdr == NULL) error("Failed to parse header. %s", args.input_fname);
    if (sam_hdr_write(args.fp_out, args.hdr)) error("Failed to write header.");
    if (args.fp_mito && bam_hdr_write(args.fp_mito, args.hdr)) error("Failed to write header.");

    // init mitochondria id
    args.mito_id = bam_name2id(args.hdr, args.mito);
    if (args.mito_id == -1) {
        warnings("Failed to find mitochondria name at the header, will NOT stat the mito ratio!");
        args.mito_id = -2; // reset to -2, in case inconsistance with unmap
    }

    // check if there is a BAM record
    if (str.s[0] != '@')
        args.preload_record = strndup(str.s, str.l);    
    free(str.s);
    return 0;
}

static void memory_release()
{
    hts_close(args.fp_out);
    ks_destroy(args.ks);
    gzclose(args.fp);
    bam_hdr_destroy(args.hdr);
    free(args.summary);    
    if (args.fp_mito) bgzf_close(args.fp_mito);
    if (args.fp_report != stdout) fclose(args.fp_report);
    if (args.enable_corr) gtf_destroy(args.G);
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

        hts_tpool *p = hts_tpool_init(nt);
        hts_tpool_process *q = hts_tpool_process_init(p, nt*2, 0);
        hts_tpool_result *r;

        for (;;) {

            struct sam_pool *b = sam_pool_read(args.ks, args.buffer_size);
            if (b == NULL) break;
            b->opts = &args;

            int block;

            do {

                block = hts_tpool_dispatch2(p, q, sam_name_parse, b, 0);

                if ((r = hts_tpool_next_result(q))) {
                    struct sam_pool *d = (struct sam_pool*)r->data;
                    write_out(d);
                }
                //hts_tpool_delete_result(r, 1);
            }
            while (block == -1);
        }
        
        hts_tpool_process_flush(q);

        while ((r = hts_tpool_next_result(q))) {
            struct sam_pool *d = (struct sam_pool *)r->data;
            write_out(d);
        }

        hts_tpool_process_destroy(q);
        hts_tpool_destroy(p);

    }

    summary_report(&args);
    
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    return 0;
}
