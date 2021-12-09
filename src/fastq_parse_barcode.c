// parse barcode sequence from reads, rename read name with barcode tags
#include "utils.h"
#include "json_config.h"
#include "ksort.h"
#include "fastq.h"
#include "kson.h"
#include "number.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "sim_search.h"

KHASH_MAP_INIT_STR(str, int)
typedef kh_str_t strhash_t;

struct bcount {
    uint64_t matched;
    uint64_t corrected;
};

struct bcode_reg {
    int rd; // read 1 or 2
    int start;
    int end;
    int dist;
    ss_t *wl;
    char **white_list; // temp allocated, will be free after initization
    int len;
    int n_wl;
};
void bcode_reg_clean(struct bcode_reg *br)
{
    if (br->wl) ss_destroy(br->wl);
}
void bcode_reg_destory(struct bcode_reg *br)
{
    bcode_reg_clean(br);
    free(br);
}
struct BRstat {
    int q30_bases;
    int bases;
    int exact_match;
    int filter;
};
struct BRstat *BRstat_alloc()
{
    struct BRstat *b = malloc(sizeof(*b));
    memset(b, 0, sizeof(*b));
    return b;
}
struct name_count_pair {
    char *name;
    uint32_t count;
};

struct fq_data {
    char *bc_str;
    int q30_bases_cell_barcode;
    int q30_bases_sample_barcode;
    int q30_bases_umi;
    int q30_bases_reads;
    int bases_cell_barcode;
    int bases_sample_barcode;
    int bases_umi;
    int bases_reads;
    int cr_exact_match;
};

#define FQ_FLAG_PASS          0
#define FQ_FLAG_BC_EXACTMATCH 1
#define FQ_FLAG_BC_FAILURE    2
#define FQ_FLAG_READ_QUAL     3
#define FQ_FLAG_SAMPLE_FAIL   4

extern int kstr_copy(kstring_t *a, kstring_t *b);

static struct args {
    const char *r1_fname;
    const char *r2_fname;
    const char *config_fname;
    const char *out1_fname;
    const char *out2_fname;
    const char *run_code;
    const char *cbdis_fname;
    const char *report_fname;
    const char *dis_fname; // barcode segment distribution

    int qual_thres;
    
    int n_thread;
    int chunk_size;
    int cell_number;

    int bgiseq_filter;
    int smart_pair;

    int dropN;
    
    // Cell barcodes found in white list, if no white list all barcodes treat as background
    struct name_count_pair *names;
    int n_name;
    int m_name;    
    strhash_t *cbhash;

    // Background
    struct name_count_pair *bgnames;
    int n_bg, m_bg;
    strhash_t *bghash;

    // file handler
    // inputs could be gzipped fastq or unzipped
    gzFile r1_fp;
    gzFile r2_fp;
    // All outputs will be unzipped for performance
    FILE *out1_fp;
    FILE *out2_fp;
    FILE *cbdis_fp;
    FILE *report_fp; // report handler
    // FILE *html_report_fp;
    FILE *barcode_dis_fp;
    
    // hold thread safe data
    // struct thread_hold **hold;
    
    // input file streaming
    struct fastq_handler *fastq;

    // stats of reads
    uint64_t raw_reads;
    uint64_t reads_pass_qc;
    uint64_t barcode_exactly_matched;
    uint64_t filtered_by_barcode;
    uint64_t filtered_by_lowqual;
    uint64_t filtered_by_sample;
    uint64_t q30_bases_cell_barcode;
    uint64_t q30_bases_sample_barcode;
    uint64_t q30_bases_umi;
    uint64_t q30_bases_reads;
    uint64_t bases_cell_barcode;
    uint64_t bases_sample_barcode;
    uint64_t bases_umi;
    uint64_t bases_reads;    
} args = {
    .r1_fname = NULL,
    .r2_fname = NULL,
    .config_fname = NULL,
    .out1_fname = NULL,
    .out2_fname = NULL,
    .run_code = NULL,
    .cbdis_fname = NULL,
    .report_fname = NULL,
    .dis_fname = NULL,
    .qual_thres = 0,
    .n_thread = 4,
    .chunk_size = 10000,
    .cell_number = 10000,
    .smart_pair = 0,
    .bgiseq_filter = 0,
    .dropN = 0,
    .names = NULL,
    .n_name = 0,
    .m_name = 0,
    .cbhash = NULL,

    .r1_fp = NULL,
    .r2_fp = NULL,
    .out1_fp = NULL,
    .out2_fp = NULL,
    .cbdis_fp = NULL,
    .report_fp = NULL,
    // .html_report_fp = NULL,
    .barcode_dis_fp = NULL,
    // .hold = NULL,
    .fastq = NULL,

    .raw_reads = 0,
    .reads_pass_qc = 0,
    .barcode_exactly_matched = 0,
    .filtered_by_barcode = 0,
    .filtered_by_lowqual = 0,
    .filtered_by_sample = 0,
    .q30_bases_cell_barcode = 0,
    .q30_bases_sample_barcode = 0,
    .q30_bases_umi = 0,
    .q30_bases_reads = 0,
    .bases_cell_barcode = 0,
    .bases_umi = 0,
    .bases_sample_barcode = 0,
    .bases_reads = 0,
};
    
static struct config {
    // consistant with white list if set
    char *cell_barcode_tag; 
    char *sample_barcode_tag;

    char *raw_cell_barcode_tag;
    char *raw_cell_barcode_qual_tag;
    int n_cell_barcode;
    struct bcode_reg *cell_barcodes;

    char *raw_sample_barcode_tag;
    char *raw_sample_barcode_qual_tag;
    int n_sample_barcode; // usually == 1, sometimes == 2,
    struct bcode_reg *sample_barcodes;

    char *umi_tag;
    char *umi_qual_tag;
    // UMI, usually random generated and located at one region
    struct bcode_reg *UMI;

    // clean read sequence
    struct bcode_reg *read_1; 
    struct bcode_reg *read_2;  // for single ends, read2 == NULL
} config = {
    .cell_barcode_tag = NULL,
    .sample_barcode_tag = NULL,
    .raw_cell_barcode_tag = NULL,
    .raw_cell_barcode_qual_tag = NULL,
    .n_cell_barcode = 0,
    .cell_barcodes = NULL,
    .raw_sample_barcode_tag = NULL,
    .raw_sample_barcode_qual_tag = NULL,
    .n_sample_barcode = 0,
    .sample_barcodes = NULL,
    .umi_tag = NULL,
    .umi_qual_tag = NULL,
    .UMI = NULL,
    .read_1 = NULL,
    .read_2 = NULL,
};
void config_destory()
{
    if (config.cell_barcode_tag) free(config.cell_barcode_tag);
    if (config.sample_barcode_tag) free(config.sample_barcode_tag);
    if (config.raw_cell_barcode_tag) free(config.raw_cell_barcode_tag);
    if (config.raw_cell_barcode_qual_tag) free(config.raw_cell_barcode_qual_tag);
    int i;
    for (i = 0; i < config.n_cell_barcode; ++i) {
        struct bcode_reg *br = &config.cell_barcodes[i];
        bcode_reg_clean(br);
    }
    free(config.cell_barcodes);
    if (config.sample_barcodes) bcode_reg_destory(config.sample_barcodes);
    if (config.UMI) bcode_reg_destory(config.UMI);
    if (config.read_1) bcode_reg_destory(config.read_1);
    if (config.read_2) bcode_reg_destory(config.read_2);
    if (config.raw_sample_barcode_tag) free(config.raw_sample_barcode_tag);
    if (config.raw_sample_barcode_qual_tag) free(config.raw_sample_barcode_qual_tag);
    if (config.umi_tag) free(config.umi_tag);
    if (config.umi_qual_tag) free(config.umi_qual_tag);
                            
}
static void config_init(const char *fn)
{
    char *config_str = json_config_open(fn);
    if (config_str == NULL) error("Empty configure file.");
    kson_t *json = kson_parse(config_str);
    free(config_str);

    const kson_node_t *root = json->root;
    if (root == NULL) error("Format error. Root node is empty.");
    int i;
    for (i = 0; i < root->n; ++i) {
        const kson_node_t *node = kson_by_index(root,i);
        if (node == NULL) continue;
        if (node->key == NULL) {
            warnings("Format error. Node key is empty. skip..");
            continue;
        }
        if (strcmp(node->key, "platform") == 0) {
            continue;
        }
        else if (strcmp(node->key, "version") == 0) {
            continue;
        }
        else if (strcmp(node->key, "cell barcode tag") == 0) {
            if (node->v.str) config.cell_barcode_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "cell barcode raw tag") == 0) {
            if (node->v.str) config.raw_cell_barcode_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "cell barcode raw qual tag") == 0) { 
            if (node->v.str) config.raw_cell_barcode_qual_tag = strdup(node->v.str);
        }        
        else if (strcmp(node->key, "cell barcode") == 0) {
            if (node->type != KSON_TYPE_BRACKET) error("Format error. \"cell barcode\":[{},{}]");
            config.n_cell_barcode = node->n;
            config.cell_barcodes = calloc(node->n,sizeof(struct bcode_reg));

            int j;
            for (j = 0; j < node->n; ++j) {
                const kson_node_t *n1 = kson_by_index(node, j);                
                if (n1 == NULL) error("cell barcode is empty.");
                if (n1->type != KSON_TYPE_BRACE) error("Format error. \"cell barcode\":[{},{}]");
                if (n1->n == 0) continue; // empty record
                struct bcode_reg *br = &config.cell_barcodes[j];
                memset(br, 0, sizeof(struct bcode_reg));
                int k;
                for (k = 0; k < n1->n; ++k) {
                    const kson_node_t *n2 = kson_by_index(n1, k);
                    if (n2 == NULL) error("cell barcode record is empty.");
                    if (strcmp(n2->key, "location") == 0) {
                        // assume format R1:1-2
                        char *p = n2->v.str;
                        if (strlen(p) < 6) error("Unknown location format, should be like \"R1:1-2\"");
                        if (p[0] == 'R' && p[2] == ':') {
                            if (p[1] == '1') br->rd = 1;
                            else if (p[1] == '2') br->rd = 2;                            
                            else  error("Unknown location format, should be like \"R[12]:1-2\"");
                        }
                        else {
                            error("Unknown location format, should be like \"R1:1-2\"");
                        }
                        p = p + 3;
                        char *s = p;
                        int c = 0;
                        for (; check_char_num(*s); s++,c++);
                        br->start = str2int_l(p, c);
                        if (*s != '-') error("Unknown location format, should be like \"R1:1-2\"");
                        p = ++s;
                        c = 0;
                        for (; check_char_num(*s); s++,c++);
                        br->end = str2int_l(p, c);
                        br->len = br->end - br->start +1;
                    }
                    else if (strcmp(n2->key, "distance") == 0) {
                        br->dist = str2int(n2->v.str);
                    }
                    else if (strcmp(n2->key, "white list") == 0) {
                        if (n2->type != KSON_TYPE_BRACKET) error("Format error. \"white list\":[]");
                        br->n_wl = n2->n;
                        br->white_list = malloc(n2->n*sizeof(char*));
                        int l;
                        for (l = 0; l < n2->n; ++l) {
                            const kson_node_t *n3 = kson_by_index(n2, l);
                            br->white_list[l] = strdup(n3->v.str);                            
                        }                                        
                    }
                    else error("Unknown key : \"%s\"", n2->key);
                }
            }
        }
        else if (strcmp(node->key, "sample barcode tag") == 0) {
            if (node->v.str) config.sample_barcode_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "sample barcode") == 0) {
            
        }
        else if (strcmp(node->key, "read 1") == 0) {
            if (node->type != KSON_TYPE_BRACE) error("Format error. \"read 1\":{}");
            config.read_1 = malloc(sizeof(struct bcode_reg));
            struct bcode_reg *br = config.read_1;
            memset(br, 0, sizeof(struct bcode_reg));
            const kson_node_t *n1 = kson_by_index(node, 0);
            if (n1 == NULL) error("read 1 location at config is empty.");
            if (strcmp(n1->key, "location") == 0) {
                // assume format R1:1-2
                char *p = n1->v.str;
                if (strlen(p) < 6) error("Unknown location format, should be like \"R1:1-2\"");
                if (p[0] == 'R' && p[2] == ':') {
                    if (p[1] == '1') br->rd = 1;
                    else if (p[1] == '2') br->rd = 2;                            
                    else  error("Unknown location format, should be like \"R[12]:1-2\"");
                }
                else {
                    error("Unknown location format, should be like \"R1:1-2\"");
                }
                p = p + 3;
                char *s = p;
                int c = 0;
                for (; check_char_num(*s); s++,c++);
                br->start = str2int_l(p, c);
                if (*s != '-') error("Unknown location format, should be like \"R1:1-2\"");
                p = ++s;
                c = 0;
                for (; check_char_num(*s); s++,c++);
                br->end = str2int_l(p, c);
                br->len = br->end - br->start +1;
            }
        }
        else if (strcmp(node->key, "read 2") == 0) {
            if (node->type != KSON_TYPE_BRACE) error("Format error. \"read 1\":{}");
            config.read_2 = malloc(sizeof(struct bcode_reg));
            struct bcode_reg *br = config.read_2;
            memset(br, 0, sizeof(struct bcode_reg));
            const kson_node_t *n1 = kson_by_index(node, 0);
            if (n1 == NULL) error("read 2 location at config is empty.");
            if (strcmp(n1->key, "location") == 0) {
                // assume format R1:1-2
                char *p = n1->v.str;
                if (strlen(p) < 6) error("Unknown location format, should be like \"R1:1-2\"");
                if (p[0] == 'R' && p[2] == ':') {
                    if (p[1] == '1') br->rd = 1;
                    else if (p[1] == '2') br->rd = 2;                            
                    else  error("Unknown location format, should be like \"R[12]:1-2\"");
                }
                else {
                    error("Unknown location format, should be like \"R1:1-2\"");
                }
                p = p + 3;
                char *s = p;
                int c = 0;
                for (; check_char_num(*s); s++,c++);
                br->start = str2int_l(p, c);
                if (*s != '-') error("Unknown location format, should be like \"R1:1-2\"");
                p = ++s;
                c = 0;
                for (; check_char_num(*s); s++,c++);
                br->end = str2int_l(p, c);
                br->len = br->end - br->start +1;
            }
        }
        else if (strcmp(node->key, "UMI tag") == 0) {
            if (node->v.str) config.umi_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "UMI qual tag") == 0) {
            if (node->v.str) config.umi_qual_tag = strdup(node->v.str);
        }
        else if (strcmp(node->key, "UMI") == 0) {
            if (node->type != KSON_TYPE_BRACE) error("Format error. \"UMI\":{}");
            config.UMI = malloc(sizeof(struct bcode_reg));
            struct bcode_reg *br = config.UMI;
            memset(br, 0, sizeof(*br));
            const kson_node_t *n1 = kson_by_index(node, 0);
            if (n1 == NULL) error("UMI location at config is empty.");
            if (strcmp(n1->key, "location") == 0) {
                // assume format R1:1-2
                char *p = n1->v.str;
                if (strlen(p) < 6) error("Unknown location format, should be like \"R1:1-2\"");
                if (p[0] == 'R' && p[2] == ':') {
                    if (p[1] == '1') br->rd = 1;
                    else if (p[1] == '2') br->rd = 2;                            
                    else  error("Unknown location format, should be like \"R[12]:1-2\"");
                }
                else {
                    error("Unknown location format, should be like \"R1:1-2\"");
                }
                p = p + 3;
                char *s = p;
                int c = 0;
                for (; check_char_num(*s); s++,c++);
                br->start = str2int_l(p, c);
                if (*s != '-') error("Unknown location format, should be like \"R1:1-2\"");
                p = ++s;
                c = 0;
                for (; check_char_num(*s); s++,c++);
                br->end = str2int_l(p, c);
                br->len = br->end - br->start +1;
            }
        }
        else {
            error("Unknown key \"%s.\"", node->key);
        }            
    }

    // init white list hash
    for (i = 0; i < config.n_cell_barcode; ++i) {
        struct bcode_reg *br = &config.cell_barcodes[i];
        if (br->n_wl == 0) continue;
        br->wl = ss_init();
        int j;
        for (j = 0; j < br->n_wl; j++) {
            int len = strlen(br->white_list[j]);
            if (len != br->len) error("Inconsistance white list length. %d vs %d, %s", br->len, len, br->white_list[j]);
            if (br->dist>3) error("Set too much distance for cell barcode, allow 3 distance at max.");
            if (br->dist > len/2) error("Allow distance greater than half of barcode! Try to reduce distance.");
            ss_push(br->wl, br->white_list[j]);
            free(br->white_list[j]);
        }
        free(br->white_list);
    }
    kson_destroy(json);
}

struct seqlite {
    kstring_t *seq;
    kstring_t *qual;
};

extern kstring_t *kstr_init();

void kstr_destory(kstring_t *str)
{
    if (str == NULL) return;
    if (str->m) free(str->s);
    free(str);
}
void seqlite_destory(struct seqlite *s)
{
    kstr_destory(s->seq);
    kstr_destory(s->qual);
    free(s);
}

struct seqlite *extract_tag(struct bseq *b, const struct bcode_reg *r, struct BRstat *stat, int *n)
{
    if (b == NULL || r == NULL) return NULL;
    
    struct seqlite *p = malloc(sizeof(*p));
    p->seq = kstr_init();
    p->qual = kstr_init();
    
    *n = 0;
    char *s = NULL;
    char *q = NULL;
    if (r->rd == 1) {
        if (r->end > b->s0.l) error("Try to select sequence out of range. [Read length: %zu, Barcode range: %d-%d. Read name: %s]", b->s0.l, r->start, r->end, b->n0.s);
        s = b->s0.s + r->start -1;
        q = b->q0.l ? b->q0.s + r->start -1 : NULL;
    }
    else {
        if (r->end > b->s1.l) error("Try to select sequence out of range. [Read length: %zu, Barcode range: %d-%d. Read name: %s]", b->s1.l, r->start, r->end, b->n0.s);
        s = b->s1.s + r->start -1;
        q = b->q1.l ? b->q1.s + r->start -1 : NULL;
    }
    int l = r->end - r->start + 1;
    kputsn(s, l, p->seq);
    kputs("", p->seq);
    if (q){
        kputsn(q, l, p->qual);
        kputs("", p->qual);
    }
    int i;
    for (i = 0; i < l; ++i) {
        if (p->qual->s[i]-33 >= 30) stat->q30_bases++;
        if (p->seq->s[i] == 'N') *n = 1;
    }
    stat->bases += l;
    return p;
}

// NULL on unfound, else on sequence
char *check_whitelist(char *s, const struct bcode_reg *r, int *exact_match)
{
    if (r->n_wl == 0) return 0;
    int len = strlen(s);
    if (len != r->len) error("Trying to check inconsistance length sequence.");
    set_hamming();
    return ss_query(r->wl, s, r->dist, exact_match);
}
static void update_rname(struct bseq *b, const char *tag, char *s){
    kstring_t str = {0,0,0};
    kputs(b->n0.s, &str);
    kputs("|||", &str);
    kputs(tag, &str);
    kputs(":Z:", &str);
    kputs(s, &str);
    kstr_copy(&b->n0, &str);
    free(str.s);
}

struct BRstat *extract_barcodes(struct bseq *b,
                                int n,
                                const struct bcode_reg *r,
                                const char *tag,
                                const char *raw_tag,
                                const char *raw_qual_tag,                                
                                const char *run_code
    )
{
    if (tag == NULL && raw_tag == NULL) return NULL;
    
    kstring_t str = {0,0,0};
    kstring_t qual = {0,0,0};
    kstring_t tag_str = {0,0,0};

    struct BRstat *stat = malloc(sizeof(struct BRstat));
    memset(stat, 0, sizeof(struct BRstat));
    stat->exact_match = 1;
    
    int dropN;
    // barcodes construct from multi segments, only all segments matched with whitelist will be export
    int i;
    for (i = 0; i < n; ++i) { 
        const struct bcode_reg *br = &r[i];
        struct seqlite *s = extract_tag(b, br, stat, &dropN);
        if (s == NULL) goto failed_check_barcode;
        char *wl = NULL;
        int exact_match = 0;
        if (br->n_wl) {
            wl = check_whitelist(s->seq->s, br, &exact_match);
            if (wl == NULL) {                
                stat->exact_match = 0;
                stat->filter = 1;
                seqlite_destory(s);
                goto failed_check_barcode;
            }
            if (exact_match == 0)  stat->exact_match = 0;
        }
        else stat->exact_match = 0;
        
        if (raw_tag) kputs(s->seq->s, &str);
        if (raw_qual_tag && s->qual->l) kputs(s->qual->s, &qual);
        
        if (tag) {            
            kputs(wl == NULL ? s->seq->s : wl, &tag_str);
        }
        seqlite_destory(s);

        if (wl) free(wl);
    }

    if (stat->filter == 1) goto failed_check_barcode;
    
    if (run_code) {
        kputc('-', &tag_str);
        kputs(run_code, &tag_str);
    }

    struct fq_data *data = (struct fq_data*)b->data;
    data->bc_str = strdup(tag_str.s);
    
    if (tag) update_rname(b, tag, tag_str.s);
    if (raw_tag) update_rname(b, raw_tag, str.s);
    if (raw_qual_tag) update_rname(b, raw_qual_tag, qual.s);

    if (str.m) free(str.s);
    if (qual.m) free(qual.s);
    if (tag_str.m) free(tag_str.s);
    return stat;
    
  failed_check_barcode:
    if (str.m) free(str.s);
    if (qual.m) free(qual.s);
    if (tag_str.m) free(tag_str.s);
    free(stat);
    return NULL;
}


struct BRstat *extract_sample_barcode_reads(struct bseq *b,
                                 int n,
                                 const struct bcode_reg *r,
                                 const char *tag,
                                 const char *raw_tag,
                                 const char *raw_qual_tag)
{
    return extract_barcodes(b, n, r, tag, raw_tag, raw_qual_tag, NULL);
}
struct BRstat *extract_cell_barcode_reads(struct bseq *b,
                                          int n,
                                          const struct bcode_reg *r,
                                          const char *tag,
                                          const char *raw_tag,
                                          const char *raw_qual_tag,                               
                                          const char *run_code)
{
    return extract_barcodes(b, n, r, tag, raw_tag, raw_qual_tag, run_code);
}

struct BRstat *extract_umi(struct bseq *b, const struct bcode_reg *r, const char *tag, const char *qual_tag)
{
    if (tag == NULL) return NULL;
    struct BRstat *stat = malloc(sizeof(*stat));
    memset(stat, 0, sizeof(*stat));

    int dropN;
    struct seqlite *s = extract_tag(b, r, stat, &dropN);

    if (args.dropN && dropN== 1) b->flag= FQ_FLAG_READ_QUAL;
    
    if (tag && s->seq->l) update_rname(b, tag, s->seq->s);
    if (qual_tag && s->qual->l) update_rname(b, qual_tag, s->qual->s);

    seqlite_destory(s);
    return stat;
}

struct BRstat *extract_reads(struct bseq *b, const struct bcode_reg *r1, const struct bcode_reg *r2)
{
    assert(r1);
    struct BRstat *stat = malloc(sizeof(struct BRstat));
    memset(stat, 0, sizeof(struct BRstat));
    int dropN;
    struct seqlite *s1 = extract_tag(b, r1, stat, &dropN);
    if (args.dropN && dropN== 1) b->flag= FQ_FLAG_READ_QUAL;

    struct seqlite *s2 = extract_tag(b, r2, stat, &dropN);
    if (args.dropN && dropN== 1) b->flag= FQ_FLAG_READ_QUAL;

    b->s0.l = 0;
    b->q0.l = 0;
    kstr_copy(&b->s0, s1->seq);
    kstr_copy(&b->q0, s1->qual);

    b->s1.l = 0;
    b->q1.l = 0;
    if (s2 && s2->seq->l) kstr_copy(&b->s1, s2->seq);
    if (s2 && s2->qual->l) kstr_copy(&b->q1, s2->qual);

    seqlite_destory(s1);
    if (s2) seqlite_destory(s2);
    
    return stat;
}

static void *run_it(void *_p)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    struct args *opts = p->opts;

    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        b->flag = FQ_FLAG_PASS;
        struct fq_data *data = malloc(sizeof(struct fq_data));
        memset(data, 0, sizeof(struct fq_data));
        b->data = data;

        if (config.sample_barcodes) {
            // sample barcode
            struct BRstat *sample_stat = extract_sample_barcode_reads(
                b,
                config.n_sample_barcode,
                config.sample_barcodes,
                config.sample_barcode_tag,
                config.raw_sample_barcode_tag,
                config.raw_sample_barcode_qual_tag);
            
            if (sample_stat == NULL) {
                b->flag = FQ_FLAG_SAMPLE_FAIL;
                continue;
            }
        
            data->q30_bases_sample_barcode = sample_stat->q30_bases;
            data->bases_sample_barcode = sample_stat->bases;
            free(sample_stat);
        }

        // UMI
        if (config.UMI) {
            struct BRstat *umi_stat = extract_umi(
                b,
                config.UMI,
                config.umi_tag,
                config.umi_qual_tag);

            if (umi_stat == NULL) error("Failed to extract UMIs.");
            
            data->q30_bases_umi = umi_stat->q30_bases;
            data->bases_umi = umi_stat->bases;
            free(umi_stat);
        }
        
        if (config.cell_barcodes) {
            // cell barcode
            struct BRstat* cell_stat = extract_cell_barcode_reads(
                b,
                config.n_cell_barcode,
                config.cell_barcodes,
                config.cell_barcode_tag,
                config.raw_cell_barcode_tag,
                config.raw_cell_barcode_qual_tag,
                opts->run_code
                );
            
            if (cell_stat == NULL) {
                b->flag = FQ_FLAG_BC_FAILURE;
                continue;
            }

            if (cell_stat->filter == 1) {
                free(cell_stat);
                b->flag = FQ_FLAG_BC_FAILURE;
                continue;
            }
            data->q30_bases_cell_barcode = cell_stat->q30_bases;
            data->bases_cell_barcode = cell_stat->bases;
            data->cr_exact_match = cell_stat->exact_match;
            if (cell_stat->exact_match) {
                b->flag = FQ_FLAG_BC_EXACTMATCH;
            }
            free(cell_stat);
        }

        if (config.read_1) {
            // clean sequence
            struct BRstat *read_stat = extract_reads(b, config.read_1, config.read_2);
            data->q30_bases_reads = read_stat->q30_bases;
            data->bases_reads = read_stat->bases;
            free(read_stat);
            if (b->flag != FQ_FLAG_PASS) continue;
            
            if (opts->bgiseq_filter) {
                if (b->q0.l) {
                    int k;
                    int bad_bases = 0;
                    for (k = 0; k < 15 && k <b->q0.l; ++k) {
                        if (b->q0.s[k]-33<10) bad_bases++;
                    }
                    if (bad_bases >2) {
                        b->flag = FQ_FLAG_READ_QUAL;
                        continue;
                    }
                    else {
                        if (b->q1.l) {
                            for (k = 0; k < 15 && k <b->q1.l; ++k) {
                                if (b->q1.s[k]-33<10) bad_bases++;
                            }
                            if (bad_bases >2) {
                                b->flag = FQ_FLAG_READ_QUAL;
                                continue;
                            }
                        }
                    }                
                }
            }

            if (opts->qual_thres > 0 && b->q0.l) { 
                int k;
                int ave = 0;
                
                for (k = 0; k < b->q0.l; k++) ave+=b->q0.s[k]-33;
                if (ave/k < opts->qual_thres) {
                    b->flag = FQ_FLAG_READ_QUAL;
                    continue;
                }

                if (b->q1.l > 0) {
                    ave = 0;
                    for (k = 0; k < b->q1.l; k++) ave+=b->q1.s[k]-33;
                    if (ave/k < opts->qual_thres) {
                        b->flag = FQ_FLAG_READ_QUAL;
                        continue;
                    }
                }
            }
        }
    }    
    return p;
}
static void write_out(void *_data)
{
    struct bseq_pool *p = (struct bseq_pool*)_data;
    struct args *opts = (struct args*)p->opts;

    FILE *fp1 = opts->out1_fp == NULL ? stdout : opts->out1_fp;
    FILE *fp2 = opts->out2_fp == NULL ? fp1 : opts->out2_fp;
    int i;
    int ret;
    // because the output queue is order, we do not consider the thread-safe of summary report
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        struct fq_data *data = (struct fq_data*)b->data;
        
        opts->raw_reads++;
        if (b->flag == FQ_FLAG_SAMPLE_FAIL) {
            opts->filtered_by_sample++;
            goto background_reads;
        }
        
        if (b->flag == FQ_FLAG_BC_FAILURE) {
            opts->filtered_by_barcode++;
            goto background_reads;
        }

        if (b->flag == FQ_FLAG_READ_QUAL) {
            opts->filtered_by_lowqual++;
            continue; // just skip ALL low quality reads
        }

        if (b->flag == FQ_FLAG_BC_EXACTMATCH) {
            opts->barcode_exactly_matched++;
            goto flag_pass;
        }

        if (b->flag == FQ_FLAG_PASS) goto flag_pass;
        
        if (0) {
          flag_pass:
            opts->reads_pass_qc++;
            fprintf(fp1, "%c%s\n%s\n", b->q0.l ? '@' : '>', b->n0.s, b->s0.s);
            if (b->q0.l) fprintf(fp1, "+\n%s\n", b->q0.s);
            if (b->s1.l > 0) {
                fprintf(fp2, "%c%s\n%s\n", b->q1.l ? '@' : '>', b->n0.s, b->s1.s);
                if (b->q1.l) fprintf(fp2, "+\n%s\n", b->q1.s);
            }

            if (data->bc_str && opts->cbhash) {
                khint_t k;
                k = kh_get(str, opts->cbhash, (char*)data->bc_str);
                if (k == kh_end(opts->cbhash)) {
                    if (opts->n_name == opts->m_name) {
                        opts->m_name += 10000;
                        opts->names = realloc(opts->names,opts->m_name *sizeof(struct name_count_pair));                        
                    }
                    struct name_count_pair *pair = &opts->names[opts->n_name];
                    pair->name = strdup(data->bc_str);
                    pair->count = 1;
                    k = kh_put(str,opts->cbhash, pair->name, &ret);
                    kh_val(opts->cbhash, k) = opts->n_name;
                    opts->n_name++;
                }
                else {
                    int id = kh_val(opts->cbhash, k);
                    struct name_count_pair *pair = &opts->names[id];
                    pair->count++;
                }
            }
        }

        if (0) {
          background_reads:                        
            if (data->bc_str && opts->bghash) {
                khint_t k;
                k = kh_get(str, opts->bghash, (char*)data->bc_str);
                if (k == kh_end(opts->bghash)) {
                    if (opts->n_bg == opts->m_bg) {
                        opts->m_bg += 10000;
                        opts->bgnames = realloc(opts->bgnames,opts->m_bg *sizeof(struct name_count_pair));                        
                    }
                    struct name_count_pair *pair = &opts->bgnames[opts->n_bg];
                    pair->name = strdup(data->bc_str);
                    pair->count = 1;
                    k = kh_put(str,opts->bghash, pair->name, &ret);
                    kh_val(opts->bghash, k) = opts->n_bg;
                    opts->n_bg++;
                }
                else {
                    int id = kh_val(opts->bghash, k);
                    struct name_count_pair *pair = &opts->bgnames[id];
                    pair->count++;
                }
            }
            
        }
        opts->q30_bases_cell_barcode += (uint64_t)data->q30_bases_cell_barcode;
        opts->q30_bases_sample_barcode += (uint64_t)data->q30_bases_sample_barcode;
        opts->q30_bases_umi += (uint64_t)data->q30_bases_umi;
        opts->q30_bases_reads += (uint64_t)data->q30_bases_reads;
        opts->bases_cell_barcode += (uint64_t)data->bases_cell_barcode;
        opts->bases_sample_barcode += (uint64_t)data->bases_sample_barcode;
        opts->bases_umi += (uint64_t)data->bases_umi;
        opts->bases_reads += (uint64_t)data->bases_reads;
        // opts->barcode_exactly_matched += data->cr_exact_match;
        if (data->bc_str) free(data->bc_str);
        free(data);
    }
    bseq_pool_destroy(p);
    fflush(fp1);
    if (fp2 != fp1) fflush(fp2);
}
static int cmpfunc (const void *a, const void *b)
{    
    return ( ((struct name_count_pair*)b)->count - ((struct name_count_pair*)a)->count );
}
void cell_barcode_count_pair_write()
{
    if (args.cbdis_fp && args.n_name && args.cbhash) {
        qsort(args.names, args.n_name, sizeof(struct name_count_pair), cmpfunc);
        int i;
        for (i = 0; i < args.n_name; ++i) {
            fprintf(args.cbdis_fp, "%s\t%d\n", args.names[i].name, args.names[i].count);
            free(args.names[i].name);
        }
        free(args.names);
        fclose(args.cbdis_fp);
    }
    if (args.cbhash) kh_destroy(str,args.cbhash);
}
void report_write()
{
    if (args.report_fp) {
        fprintf(args.report_fp, "Number of Fragments,%"PRIu64"\n", args.raw_reads);
        fprintf(args.report_fp, "Fragments pass QC,%"PRIu64"\n", args.reads_pass_qc);
        fprintf(args.report_fp, "Fragments with Exactly Matched Barcodes,%"PRIu64"\n", args.barcode_exactly_matched);
        fprintf(args.report_fp, "Fragments with Failed Barcodes,%"PRIu64"\n", args.filtered_by_barcode);
        fprintf(args.report_fp, "Fragments Filtered on Low Quality,%"PRIu64"\n", args.filtered_by_lowqual);
        fprintf(args.report_fp, "Fragments Filtered on Unknown Sample Barcodes,%"PRIu64"\n", args.filtered_by_sample);
        fprintf(args.report_fp, "Q30 bases in Cell Barcode,%.1f%%\n",args.bases_cell_barcode == 0 ? 0 : (float)args.q30_bases_cell_barcode/(args.bases_cell_barcode+1)*100);
        fprintf(args.report_fp, "Q30 bases in Sample Barcode,%.1f%%\n", args.bases_sample_barcode == 0 ? 0 : (float)args.q30_bases_sample_barcode/(args.bases_sample_barcode+1)*100);
        fprintf(args.report_fp, "Q30 bases in UMI,%.1f%%\n", args.bases_umi == 0 ? 0 : (float)args.q30_bases_umi/(args.bases_umi+1)*100);
        fprintf(args.report_fp, "Q30 bases in Reads,%.1f%%\n", (float)args.q30_bases_reads/(args.bases_reads+1)*100);
        
        if (args.report_fp != stderr) fclose(args.report_fp);
    }
}
void full_details()
{
    /*
    if (args.barcode_dis_fp) {
        LOG_print("Cell barcodes summary.");
        int i, j, k;
      
        for (i = 1; i < args.n_thread; ++i) {
            for (j = 0; j < args.hold[i]->n; ++j) {
                struct segment *s0 = args.hold[0]->seg[j];
                struct segment *s1 = args.hold[i]->seg[j];
                for (k = 0; k < s1->n; ++k) {
                    s0->counts[k].matched += s1->counts[k].matched;
                    s0->counts[k].corrected += s1->counts[k].corrected;
                }
            }
        }
      
        for (i = 0; i < config.n_cell_barcode; ++i) {
            struct bcode_reg *br = &config.cell_barcodes[i];
            struct segment *s0 = args.hold[0]->seg[i];
            fprintf(args.barcode_dis_fp, "# Read %d, %d-%d\n", br->rd, br->start, br->end);
            for (j = 0; j < br->n_wl; ++j)
                fprintf(args.barcode_dis_fp, "%s\t%"PRIu64"\t%"PRIu64"\n", br->white_list[j], s0->counts[j].matched, s0->counts[j].corrected);
        }
    }
    */

}
static void memory_release()
{
    if (args.r1_fp) gzclose(args.r1_fp);
    if (args.r2_fp) gzclose(args.r2_fp);
    if (args.out1_fp) fclose(args.out1_fp);
    if (args.out2_fp) fclose(args.out2_fp);
    if (args.barcode_dis_fp) fclose(args.barcode_dis_fp);
    fastq_handler_destory(args.fastq);
}
static int parse_args(int argc, char **argv)
{
    if ( argc == 1 ) return 1;

    int i;
    const char *thread = NULL;
    const char *chunk_size = NULL;    
    const char *qual_thres = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;

        if (strcmp(a, "-1") == 0) var = &args.out1_fname;
        else if (strcmp(a, "-2") == 0) var = &args.out2_fname;
        else if (strcmp(a, "-config") == 0) var = &args.config_fname;
        else if (strcmp(a, "-cbdis") == 0) var = &args.cbdis_fname;
        else if (strcmp(a, "-t") == 0) var = &thread; // skip
        else if (strcmp(a, "-r") == 0) var = &chunk_size; // skip
        else if (strcmp(a, "-run") == 0) var = &args.run_code;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-dis") == 0) var = &args.dis_fname;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;       
        else if (strcmp(a, "-f") == 0) {
            args.bgiseq_filter = 1;
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

    if (args.config_fname == NULL) error("Option -config is required.");
    config_init(args.config_fname);
    LOG_print("Configure file inited.");
    
    if (thread) args.n_thread = str2int((char*)thread);
    //if (chunk_size) args.chunk_size = str2int((char*)chunk_size);
    // assert(args.n_thread >= 1 && args.chunk_size >= 1);
    if (qual_thres) {
        args.qual_thres = str2int((char*)qual_thres);
        LOG_print("Average quality below %d will be drop.", args.qual_thres);
    }

    if (args.r1_fname == NULL && (!isatty(fileno(stdin)))) args.r1_fname = "-";
    if (args.r1_fname == NULL) error("Fastq file(s) must be set.");
        
    if (args.report_fname) {
        args.report_fp = fopen(args.report_fname, "w");
        CHECK_EMPTY(args.report_fp, "%s : %s.", args.report_fname, strerror(errno));
    }
    else args.report_fp = stderr;

    if (args.out1_fname) {
        args.out1_fp = fopen(args.out1_fname, "w");
        if (args.out1_fp == NULL) error("%s: %s.", args.out1_fname, strerror(errno));
        if (args.out2_fname) {
            args.out2_fp = fopen(args.out2_fname,"w");
            if (args.out2_fp == NULL) error("%s: %s.", args.out2_fname, strerror(errno));
        }
    }
    
    args.fastq = fastq_handler_init(args.r1_fname, args.r2_fname, args.smart_pair, args.chunk_size);
    if (args.fastq == NULL) error("Failed to init input fastq.");

    if (args.cbdis_fname) {
        args.cbdis_fp = fopen(args.cbdis_fname, "w");
        if (args.cbdis_fp == NULL) error("%s : %s.", args.cbdis_fname, strerror(errno));
        args.cbhash = kh_init(str);
        args.bghash = kh_init(str);
    } 
        
    // if (args.run_code == NULL) args.run_code = strdup("1");

    if (args.dis_fname != NULL) {
        args.barcode_dis_fp = fopen(args.dis_fname, "w");
        CHECK_EMPTY(args.barcode_dis_fp, "%s : %s.", args.dis_fname, strerror(errno));
    }
    
    return 0;
}

extern int fastq_parse_usage();
static void *process(void *shared, int step, void *_data)
{
    struct args *aux = (struct args *)shared;
    struct bseq_pool *data = (struct bseq_pool*)_data;

    if (step == 0) {
        struct bseq_pool *b = fastq_read(aux->fastq, &args);
        return b;
    }
    else if (step == 1) {
        return run_it(data);        
    }
    else if (step == 2) {
        write_out(data);
        return 0;
    }
    return 0;
}
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

int fastq_prase_barcodes(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return fastq_parse_usage();

    if (args.n_thread == 1) {
        for (;;) {
            struct bseq_pool *b = fastq_read(args.fastq, &args);
            if (b == NULL) break;
            b = run_it(b);
            write_out(b);
        }
    }

    kt_pipeline(args.n_thread, process, &args, 3);
    
    cell_barcode_count_pair_write();

    report_write();
    
    // full_details();
    
    memory_release();

    config_destory();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());    
    return 0;
}
