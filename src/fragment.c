#include "utils.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "dict.h"
#include "number.h"
#include "region_index.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *barcode_list;
    const char *sites_fname;
    const char *tag;
    int isize;
    int file_th;
    int qual_thres;
    struct dict *cells;
    int disable_offset;
    int reverse_offset;
    int forward_offset;
} args = {
    .input_fname    = NULL,
    .output_fname   = NULL,
    .barcode_list   = NULL,
    .sites_fname    = NULL,
    .tag            = NULL,
    .isize          = 2000,
    .file_th        = 4,
    .qual_thres     = 20,
    .cells          = NULL,
    .disable_offset = 0,
    .reverse_offset = -5,
    .forward_offset = 4,
};

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *file_th = NULL;
    const char *isize   = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-list") == 0) var = &args.barcode_list;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-isize") == 0) var = &isize;
        else if (strcmp(a, "-disable-offset") == 0) {
            args.disable_offset = 1;
            continue;
        }
        if (var != 0) {
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        
        error("Unknown argument : %s", a);
    }

    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (args.output_fname == NULL) error("No fragment file set.");
    if (args.tag == NULL) {
        warnings("-tag is not set. Trying to use \"CB\" tag for cell barcode.");
        args.tag = "CB";
    }
    if (strlen(args.tag) != 2) error("Bad format of tag, %s", args.tag);
    if (file_th) args.file_th = str2int(file_th);
    if (isize) args.isize = str2int(isize);
    if (args.isize < 0) args.isize = 2000;
    
    args.cells = dict_init();
    dict_set_value(args.cells);

    if (args.barcode_list)
        dict_read(args.cells, args.barcode_list);

    if (args.disable_offset) {
        args.reverse_offset = 0;
        args.forward_offset = 0;
    }
    return 0;
}  

struct frag {
    int tid;
    int start;
    int end;
    int dup;
    int idx;
    struct frag *next;
};

struct frag_pool {
    struct frag *head;
    struct frag *tail;
    int n;         // cached nodes
    int cut_sites; // count of all nodes
    struct region_index *idx;
};
void export_sites_stat(struct dict *d, const char *sites_fname)
{
    if (sites_fname == NULL) return;
    
    FILE *fp = fopen(sites_fname, "w");
    if (fp == NULL) error("%s: %s", sites_fname, strerror(errno));
    int i;
    for (i = 0; i < dict_size(d); ++i) {
        struct frag_pool *p = dict_query_value(d, i);
        if (p == NULL) continue;
        if (p->n == 0) continue;
        fprintf(fp, "%s\t%d\n", dict_name(d, i), p->cut_sites);
    }
    fclose(fp);
}
struct frag_pool *fragment_pool_create()
{
    struct frag_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    p->idx = region_index_create();
    return p;
}
static int cmpfunc(const void *_a, const void *_b)
{
    struct frag *a = *(struct frag**)_a;
    struct frag *b = *(struct frag**)_b;
    return a->start - b->start == 0 ? a->end - b->end : a->start - b->start;
}

void fragment_flush_cache(struct dict *d, BGZF *out, bam_hdr_t *hdr)
{
    assert(out);
    
    int i;
    int sizes = 0;
    void **a = NULL;
    for (i = 0; i < dict_size(d); ++i) {
        struct frag_pool *p = dict_query_value(d, i);
        if (p == NULL) continue;
        debug_print("%d", p->n);
        if (p->n == 0) continue;
        int j = sizes;
        sizes += p->n;
        a = realloc(a, sizeof(void**)*sizes);
        struct frag *p0 = p->head;
        while (p0) {
            a[j++] = p0;
            p0 = p0->next;
        }
        assert(j == sizes);
    }
    
    qsort(a, sizes, sizeof(struct frag*), cmpfunc);

    // write to disk
    kstring_t str = {0,0,0};
    for (i = 0; i < sizes; ++i) {
        str.l = 0;
        struct frag *f = a[i];
        //debug_print("%p",f);
        debug_print("%d", f->tid);
        kputs(hdr->target_name[f->tid], &str); kputc('\t', &str);        
        kputw(f->start, &str); kputc('\t', &str);
        kputw(f->end, &str); kputc('\t', &str);
        kputs(dict_name(d, f->idx), &str); kputc('\t', &str);
        kputw(f->dup, &str);kputs("\n", &str);
        if (bgzf_write(out, str.s, str.l) < 0) error("Failed to write file.");
    }
    free(str.s);

    // release cached
    for (i = 0; i < dict_size(d); ++i) {
        struct frag_pool *p = dict_query_value(d, i);
        if (p == NULL) continue;
        if (p->n == 0) continue;
        struct frag *p0 = p->head;
        while (p0) {
            p->n--;
            p->head = p0->next;
            free(p0);
            p0 = p->head;
        }
        p->head = p->tail = NULL;
        assert(p->n == 0);
        region_index_destroy(p->idx);
    }
    
}
void fragment_close(struct dict *d)
{
    int i;
    for (i = 0; i < dict_size(d); ++i) {
        struct frag_pool *p = dict_query_value(d, i);        
        if (p == NULL) continue;
        if (p->n == 0) continue;
        error("%s is still reachable.", dict_name(d, i));
        free(p);
    }
    dict_destroy(d);
}
void fragment_pool_push0(struct frag_pool *p, int start, int end, int tid, int idx)
{
    if (p->idx == NULL)
        p->idx = region_index_create(); // reset

    struct region_itr *itr = region_query(p->idx, start, end);
    
    int i;
    for (i = 0; itr && i < itr->n; ++i) {
        struct frag *f0 = (struct frag*)itr->rets[i];
        if (f0->start == start && f0->end == end) {
            f0->dup++;
            return; // duplication
        }
    }

    struct frag *f = malloc(sizeof(*f));
    memset(f, 0, sizeof(*f));
    f->start = start;
    f->end = end;
    f->tid = tid;
    f->idx = idx;
    if (p->head && p->tail) {
        p->tail->next = f;
        p->tail = f;
    }
    else {
        p->head = p->tail = f;
    }
    index_bin_push(p->idx, start, end, f);
    p->n++;
    p->cut_sites++;
}
// 1 on failure, 0 on success
static int fragment_pool_push(struct dict *cells, bam1_t *b, const char *CB, int wl)
{
    uint8_t *val = bam_aux_get(b, CB);
    if (!val) return 1;
    
    int ret = dict_query(cells, (char*)(val+1));
    if (wl == 1 && ret == -1) return 1;
    if (ret == -1) ret = dict_push(cells, (char*)(val+1));
    
    int start, end;
    if (b->core.isize > 0) {
        start = b->core.pos+1;
        end = start + b->core.isize;
    }
    else {
        end = b->core.pos+1;
        start = end + b->core.isize;
    }

    // Tn5 offset
    if (b->core.flag & BAM_FREVERSE) {
        start = start + args.reverse_offset;
        end   = end + args.reverse_offset;
    }
    else {
        start = start + args.forward_offset;
        end   = end + args.forward_offset;
    }

    struct frag_pool *p = dict_query_value(cells, ret);
    if (p == 0) {
        p = fragment_pool_create();
        dict_assign_value(cells, ret, p);
    }
    fragment_pool_push0(p, start, end, b->core.tid, ret);

    return 0;
}

extern int fragment_usage();

int bam2frag(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return fragment_usage();

    htsFile *fp  = hts_open(args.input_fname, "r");
    if (fp == NULL)
        error("%s : %s.", args.input_fname, strerror(errno));
    
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");

    hts_set_threads(fp, args.file_th);

    BGZF *fp_out = bgzf_open(args.output_fname, "w");
    if (fp_out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    
    bgzf_mt(fp_out, args.file_th, 256);

    bam1_t *b = bam_init1();
    int last_id = -1;
    int ret;
    while ((ret = sam_read1(fp, hdr, b)) >=0) {
        if (args.qual_thres > 0 && b->core.qual < args.qual_thres) continue;
        if (b->core.tid < 0) continue;
        if (b->core.flag & BAM_FREAD2) continue; 
        if (last_id == -1) last_id = b->core.tid;
        
        // output buffered Records
        if (last_id != b->core.tid) {
            fragment_flush_cache(args.cells, fp_out, hdr);
            last_id = b->core.tid;
        }

        fragment_pool_push(args.cells, b, args.tag, args.barcode_list ? 1 : 0);
    } 
    
    fragment_flush_cache(args.cells, fp_out, hdr);

    export_sites_stat(args.cells,args.sites_fname);

    fragment_close(args.cells);    

    bgzf_close(fp_out);     
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(fp);

    LOG_print("Fragment file created. Indexing..");
    // build index
    const tbx_conf_t tbx_conf_bed = { TBX_UCSC, 1, 2, 3, '#', 0 };
    
    if (tbx_index_build(args.output_fname, 0, &tbx_conf_bed))
        warnings("Failed to build index file of %s.", args.output_fname);

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
