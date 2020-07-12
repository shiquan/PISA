#include "utils.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "bam_pool.h"
#include "htslib/kstring.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include "thread.h"
#include "dict.h"
#include "htslib/khash.h"

// maximal hamming distance of two similar UMI
#define UMI_E 1

extern char *compactDNA(const char *a, int l);
extern int   compDNA_hamming_distance(const char *a, const char *b);
extern char *compDNA_decode(const char *a);

typedef struct umi_count {
    int count; 
    int index; // index refer to umi dict
    int umi_index; // index to similiar UMIs group
    int primary;
    int filter;
} uc_t;

KHASH_MAP_INIT_STR(bc, uc_t)

// To correct UMI needs first cache all UMIs and grouped by predefined barcodes, such as cell and gene.
// For example, if set -tags-block to "CB,GN", this will group reads from the same cell barcode and gene
// annotation first, reads in one group will future be grouped by UMIs, if two subgroups have similar
// UMIs (==1 hamming distance by default), the UMI of less supported is corrected to the UMI with higher support.

// This cache structure is a dictionary structure, use cell barcode as the key and point to "struct bc_corr".
// The "struct bc_corr" is an iterative structure, point to itself if more than 1 tag defined by -tags-block.
// cell::bc_corr::tag::tag::UMIs
struct tag_val {
    struct dict *bc;
    kh_bc_t     *val;
};

struct bc_corr {
    struct dict    *umi_val; // UMIs
    struct tag_val *val;
};

void bc_corr_destroy1(struct dict *bc)
{
    int i;
    for (i = 0; i < dict_size(bc); ++i) {
        struct tag_val *v = dict_query_value(bc, i);
        if (v->bc != NULL)
            bc_corr_destroy1(v->bc);
        else 
            kh_destroy(bc, v->val);
    }
    dict_destroy(bc);
}
void bc_corr_destroy(struct dict *C)
{
    int i;
    for (i = 0; i < dict_size(C); ++i) {
        struct bc_corr *bc = dict_query_value(C, i);
        if (bc->val->bc) 
            bc_corr_destroy1(bc->val->bc);
        else 
            kh_destroy(bc,bc->val->val);
        free(bc->val);
        dict_destroy(bc->umi_val);
        free(bc);
    }
    dict_destroy(C);
}

// this function return a array of sam attributions, values of this array are actually point to SAM::data, so don't free them
char **sam_tag_values(bam1_t *b, int n, const char **blocks)
{
    char **v = malloc(n*sizeof(void*));
    int i;
    for (i = 0; i <n; ++i) {
        char *v0 = (char*)bam_aux_get(b, blocks[i]);
        if (!v0) {
            free(v);
            return NULL;
        }
        v[i] = v0+1; // skip the type character
    }
    return v;
}

static struct args {
    const char  * input_fname;
    const char  * output_fname;

    const char  * tag;
    int           n_block;
    char        **blocks;
    
    const char  * new_tag;

    int           cr_method; // correct umi like cellrange
    
    htsFile     * in;
    htsFile     * out;

    bam_hdr_t   * hdr;
    
    int           n_thread;
    int           file_th;
    
    int           chunk_size;
    struct dict * Cindex;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .tag          = NULL,
    .new_tag      = NULL,
    .n_block      = 0,
    .blocks       = NULL,
    .cr_method    = 0,
    .in           = NULL,
    .out          = NULL,
    .hdr          = NULL,
    .n_thread     = 5,
    .file_th      = 4,

    .chunk_size   = 1000000, //1M
    .Cindex       = NULL,
};

static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.in);
    sam_close(args.out);
    int i;
    for (i = 0; i < args.n_block; ++i) free(args.blocks[i]);
    free(args.blocks);
    bc_corr_destroy(args.Cindex);
}
static void umi_idx_refresh(kh_bc_t *val, int old_idx, int new_idx)
{
    khiter_t k;
    for (k = kh_begin(val); k != kh_end(val); ++k) {
        if (!kh_exist(val, k)) continue;       
        struct umi_count *cnt = &kh_val(val, k);
        if (cnt->umi_index == old_idx) cnt->umi_index = new_idx;
    }
}
static void umi_idx_correct(kh_bc_t *val)
{
    khiter_t k;
    for (k = kh_begin(val); k != kh_end(val); ++k) {
        if (!kh_exist(val, k)) continue;       
        struct umi_count *cnt = &kh_val(val, k);
        if (cnt->umi_index <= 0) continue; // corrected or no need to correct

        int best_umi = cnt->index;
        int best_cnt = cnt->count;

        // check the best UMI
        khiter_t k1;
        for (k1 = k+1; k1 != kh_end(val); ++k1) {
            if (!kh_exist(val, k1)) continue;
            struct umi_count *cnt1 = &kh_val(val, k1);
            if (cnt1->umi_index != cnt->umi_index) continue; // not in same UMI group
            if (best_cnt < cnt1->count) {
                best_umi = cnt1->index;
                best_cnt = cnt1->count;
            }
            
        }

        // update UMI index
        for (k1 = k+1; k1 != kh_end(val); ++k1) {
            if (!kh_exist(val, k1)) continue;
            struct umi_count *cnt1 = &kh_val(val, k1);
            if (cnt1->umi_index != cnt->umi_index) continue;
            if (cnt1->index != best_umi) cnt1->primary = 0; // not primary UMI
            cnt1->index = best_umi;
            cnt1->umi_index = -1; // checked
        }
        if (cnt->index != best_umi) cnt->primary = 0;
        cnt->index  = best_umi;
        cnt->umi_index = -1; // checked

    }
}

void build_index_core(struct tag_val *tag_val, struct dict *umi_val)
{
    if (tag_val->bc != NULL) { //iterate next tag
        int i;
        for (i = 0; i < dict_size(tag_val->bc); ++i) {
            char *name = dict_name(tag_val->bc, i);
            int id = dict_query(tag_val->bc, name);
            struct tag_val *v = dict_query_value(tag_val->bc, id);
            build_index_core(v, umi_val);
        }
        return;
    }

    // UMI correction, O(n^2)
    kh_bc_t *val = tag_val->val;
    khiter_t k;
    
    int umi_idx = 1;

    for (k = kh_begin(val); k != kh_end(val); ++k) {
        if (!kh_exist(val, k)) continue;
        
        struct umi_count *cnt = &kh_val(val, k);
        const char *a = kh_key(val, k);
        
        khiter_t k1;
        for (k1 = k+1; k1 != kh_end(val); ++k1) {
            if (!kh_exist(val, k1)) continue;            
            struct umi_count *cnt1 = &kh_val(val,k1);
            if (cnt1->umi_index == cnt->umi_index && cnt->umi_index > 0) continue; // umi_idx already updated
            if (cnt1->umi_index == umi_idx) continue; // umi_idx already updated and refreshed            
            
            // check similarity of two UMIs
            const char *b = kh_key(val, k1);
            int e = compDNA_hamming_distance(a, b);
            if (e > UMI_E) continue;
            
            if (cnt->umi_index == 0) cnt->umi_index = umi_idx; // init the UMI index; if no similar UMI, umi_index == 0
            
            if (cnt1->umi_index == 0) cnt1->umi_index = cnt->umi_index;
            if (cnt1->umi_index != umi_idx) // in case already point to another group
                umi_idx_refresh(val, cnt1->umi_index, cnt->umi_index); // update idx to new index
        }
        
        umi_idx++; // increase the index flag
    }

    umi_idx_correct(val);
}
// for one cell, reads with the same umi but map to more than one gene, only keep gene with higher read support.
// in case of a tie for maximal read support, all reads are discarded.
void filter_umi_gene(struct bc_corr *bc)
{
    assert(bc->val->bc); // check if gene exists
    dict_set_value(bc->umi_val);
    struct dict *gene = bc->val->bc;
    int i;
    for (i = 0; i < dict_size(gene); ++i) {
        struct tag_val *v = dict_query_value(gene, i);
        assert(v->val);
        
        khiter_t k;
        for (k = kh_begin(v->val); k != kh_end(v->val); ++k) { // for each UMI in one gene
            if (!kh_exist(v->val, k)) continue;
            
            struct umi_count *cnt = &kh_val(v->val,k);
            if (cnt->primary == 0) continue; // skip if not primary UMI
            
            struct umi_count *cnt0 = dict_query_value(bc->umi_val, cnt->index);
            if (cnt0 == NULL) 
                dict_assign_value(bc->umi_val, cnt->index, cnt);
            else { // already present
                
                // compare the read count of two group
                if (cnt0->count < cnt->count) {
                    cnt0->filter = 1;
                    dict_assign_value(bc->umi_val, cnt->index, cnt); // keep gene UMI with high frequency
                }
                else if (cnt0->count == cnt->count) {
                    cnt0->filter = 1; // ambigous
                    cnt->filter = 1;
                }
                else {
                    cnt->filter = 1;
                }
            }
        }
    }
}
void build_index1(struct bc_corr *bc)
{
    build_index_core(bc->val, bc->umi_val);
    if (args.cr_method)
        filter_umi_gene(bc);
}

kh_bc_t *select_umi_hash(struct dict *Cindex, int n, const char **tags)
{
    int cell_idx = dict_query(Cindex, tags[0]);
    if (cell_idx < 0) cell_idx = dict_push(Cindex, tags[0]);
    struct bc_corr *bc = dict_query_value(Cindex, cell_idx);
    if (bc == NULL) {
        bc = malloc(sizeof(struct bc_corr));
        bc->val = malloc(sizeof(struct tag_val));
        bc->umi_val = dict_init();
        
        bc->val->val = NULL;
        bc->val->bc = NULL;
        dict_assign_value(Cindex, cell_idx, bc);
    }
    
    struct tag_val *v = bc->val;
    int i;    
    for (i = 1; i < n; ++i) {
        if (v->bc == NULL) {
            v->bc = dict_init();
            dict_set_value(v->bc);
        }
        int ret = dict_query(v->bc, tags[i]);
        if (ret < 0) ret = dict_push(v->bc, tags[i]);
        
        struct tag_val *v0 = dict_query_value(v->bc, ret);
        if (v0 == NULL) {
            v0 = malloc(sizeof(*v0));
            v0->bc = NULL;
            v0->val = NULL;
            dict_assign_value(v->bc, ret, v0);            
        }
        v = v0; // update point to next level
    }

    if (v->val == NULL) v->val = kh_init(bc);
    return v->val;
}
void bc_push(struct dict *bc, int cr_method, int n_tag, const char **tags, const char *umi_tag, bam1_t *b)
{
    char **v = sam_tag_values(b, n_tag, tags);
    if (v == NULL) return;
    char *umi = (char*)bam_aux_get(b, umi_tag);
    if (!umi) {
        free(v);
        return;
    }

    umi = umi+1; // emit type
    
    kh_bc_t *uhash = select_umi_hash(bc, n_tag, (const char**)v); // select_umi_hash will auto init empty cell
    char *comp = compactDNA(umi, strlen(umi));
    
    struct bc_corr *bc0 = dict_query_value2(bc, v[0]);
    free(v);
    int id = dict_push(bc0->umi_val, comp);
    free(comp);
    char *cc = dict_name(bc0->umi_val, id); // use cached string

    khiter_t k;    
    k = kh_get(bc, uhash, cc);
    if (k == kh_end(uhash)) {
        int ret;
        k = kh_put(bc, uhash, cc, &ret);
        struct umi_count *uc = &kh_val(uhash, k);
        uc->count     = 1;
        uc->index     = id; // index of umi dict
        uc->umi_index = 0;  // init state, do NOT change it
        uc->filter    = 0;  // init state, do NOT change it
        uc->primary   = 1;  // any UMI is primary before check
    }
    else {
        struct umi_count *uc = &kh_val(uhash, k);
        uc->count++;
    }
}

struct dict *build_index(const char *fn, int cr_method, int n_tag, const char **tags, const char *umi_tag)
{
    LOG_print("Building index ..");
    double t_real;
    t_real = realtime();

    struct dict *cell_bc = dict_init();
    dict_set_value(cell_bc);
    
    htsFile *fp = hts_open(fn, "r");
    
    CHECK_EMPTY(fp, "%s : %s.", fn, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam) error("Unsupported input format, only support BAM/SAM/CRAM format.");
    
    bam_hdr_t *hdr = sam_hdr_read(fp);
    hts_set_threads(fp, args.file_th);
    bam1_t *b = bam_init1();
    int i;
    for (;;) {
        if (sam_read1(fp, hdr, b) < 0) break;
        bam1_core_t *c = &b->core;
        if (c->flag & BAM_FQCFAIL ||
            c->flag & BAM_FSECONDARY ||
            c->flag & BAM_FSUPPLEMENTARY ||
            c->flag & BAM_FUNMAP ||
            c->flag & BAM_FDUP) continue;
        bc_push(cell_bc, cr_method, n_tag, tags, umi_tag, b);
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);

    //struct tpool *tp = tpool_init(args.n_thread, args.n_thread*2, 1);

    for (i = 0; i < dict_size(cell_bc); ++i) {
        struct bc_corr *bc0 = dict_query_value(cell_bc, i);
        build_index1(bc0);
        //tpool_add_work(tp, build_index1, bc0);
    }
    //tpool_destroy(tp);
    LOG_print("Build time : %.3f sec", realtime() - t_real);
    return cell_bc;
}

static int parse_args(int argc, char **argv)
{
    if (argc ==1) return 1;
    
    const char *block_tags = NULL;
    const char *file_th = NULL;
    const char *thread = NULL;
    
    int i;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-tags-block") == 0) var = &block_tags;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-new-tag") == 0) var = &args.new_tag;
        else if (strcmp(a, "-cr") == 0) {
            args.cr_method = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1] != '\0') error("Unknown option, %s", a);
        
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument : %s",a);
    }

    CHECK_EMPTY(args.output_fname, "-o need to be set.");
    CHECK_EMPTY(args.input_fname, "No input bam.");
    CHECK_EMPTY(args.tag, "-tag need to be set.");
    CHECK_EMPTY(block_tags, "-tags-block need to be set.");
    
    kstring_t str = {0,0,0};
    kputs(block_tags, &str);
    int *s = ksplit(&str, ',', &args.n_block);
    assert(args.n_block >0);
    args.blocks = malloc(args.n_block*sizeof(char*));
    for (i = 0; i < args.n_block; ++i)
        args.blocks[i] = strdup(str.s+s[i]);
    free(s);
    free(str.s);

    if (file_th) args.file_th = str2int((char*)file_th);
    if (thread) args.n_thread = str2int((char*)thread);

    args.Cindex = build_index(args.input_fname, args.cr_method, args.n_block, (const char **)args.blocks, args.tag);

    if (args.Cindex == NULL) error("Failed to index.");
    
    return 0;
}
char *select_umi(struct dict *Cindex, int n,const char **tags, char *umi)
{
    int l = strlen(umi);
    char *cu = compactDNA(umi, l);
    assert(cu);
    kh_bc_t *v = select_umi_hash(Cindex, n, tags);
    khiter_t k;
    k = kh_get(bc, v, cu);
    assert(k != kh_end(v));
    free(cu);
    struct umi_count *cnt = &kh_val(v,k);
    if (cnt->filter) return NULL;
    
    int cell_idx = dict_query(Cindex, tags[0]);
    struct bc_corr *bc = dict_query_value(Cindex, cell_idx);
    const char *r = dict_name(bc->umi_val, cnt->index);
    return compDNA_decode(r);
}
int update_new_tag(struct dict *Cindex, int n_block, const char **blocks, const char *old_tag, const char *new_tag, bam1_t *b)
{
    char *umi = (char *)bam_aux_get(b, old_tag);
    assert(umi);
    
    char **tag_vals = sam_tag_values(b, n_block, blocks);
    if (tag_vals == NULL) return 0;
    
    char *new_umi = select_umi(Cindex, n_block, (const char **)tag_vals, umi+1);
   
    free(tag_vals);
    
    if (!new_umi) return 0;
    
    if (new_tag)
        bam_aux_append(b, args.new_tag, 'Z', strlen(new_umi)+1, (uint8_t*)new_umi);
    else
        memcpy(umi, new_umi, strlen(new_umi)); // since it is equal length, just reset the memory..
    
    free(new_umi);
    return 1;
}

static void *run_it(void *data)
{
    struct bam_pool *p = (struct bam_pool*)data;
    int i;
    int c = 0;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = &p->bam[i];
        c += update_new_tag(args.Cindex, args.n_block, (const char**)args.blocks, args.tag, args.new_tag, b);
    }
    
    return p;
}
static void write_out(struct bam_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i)        
        if (sam_write1(args.out, args.hdr, &p->bam[i]) == -1) error("Failed to write SAM.");
    bam_pool_destory(p);
}

extern int bam_corr_usage();

int bam_corr_umi(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return bam_corr_usage();

    args.in  = hts_open(args.input_fname, "r");
    args.hdr = sam_hdr_read(args.in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
    
    hts_set_threads(args.out, args.file_th); // write file in multi-threads

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    for (;;) {
        struct bam_pool *b = bam_pool_create();
        bam_read_pool(b, args.in, args.hdr, args.chunk_size);
            
        if (b == NULL) break;
        if (b->n == 0) { free(b->bam); free(b); break; }
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
                write_out(d);   
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);        
    }
    
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    
    memory_release();    
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    //LOG_print("%d records updated.", args.update_count);
    return 0;
}
    
#ifdef CORR_UMI
int main(int argc, char **argv)
{
    return bam_corr_umi(argc, argv);
}

#endif
