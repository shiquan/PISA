// parse barcode sequence from reads, rename read name with barcode tags
#include "utils.h"
#include "thread_pool.h"
#include "json_config.h"
#include "ksort.h"
#include "fastq.h"
#include "kson.h"
#include "number.h"
#include <pthread.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"

KHASH_MAP_INIT_STR(str, int)
typedef kh_str_t strhash_t;

const char *program_name = "CellBarcodeParser";

static int usage()
{
    fprintf(stderr, "* Parse cell barcode and UMI string from raw FASTQ.\n");
    fprintf(stderr, " %s [options] read_1.fq.gz read_2.fq.gz\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -1      [fastq]    Read 1 output.\n");
    fprintf(stderr, "  -2      [fastq]    Read 2 output.\n");
    fprintf(stderr, "  -config [txt]      Configure file in JSON format. Required.\n");
    fprintf(stderr, "  -run    [string]   Run code, used for different library.\n");
    fprintf(stderr, "  -cbdis  [txt]      Cell barcode sequence and count pairs.\n");
    fprintf(stderr, "  -t      [INT]      Thread.\n");
    fprintf(stderr, "  -r      [INT]      Records per chunk. [10000]\n");
    fprintf(stderr, "  -report [txt]      Summary report.\n");
    fprintf(stderr, "  -dis    [txt]      Barcode distribution count.\n");
    fprintf(stderr, "  -f                 Filter reads based on BGISEQ standard. Two bases quality < q10 at first 15.\n");
    fprintf(stderr, "  -q      [INT]      Drop this read if average sequencing quality below this value.\n");
    fprintf(stderr, "\n");
    return 1;
}

struct bcount {
    uint64_t matched;
    uint64_t corrected;
};

struct segment {
    int n; // white list barcodes per segment
    struct bcount *counts;
};
struct segment *segment_alloc(int n)
{
    struct segment *s = malloc(sizeof(*s));
    s->n = n;
    s->counts = malloc(sizeof(struct bcount)*n);
    memset(s->counts, 0, sizeof(struct bcount)*n);
    return s;
}
struct thread_hold {
    int n; // segments
    struct segment **seg;
};
struct thread_hold *thread_hold_fork(struct thread_hold *h)
{
    struct thread_hold *f = malloc(sizeof(*f));
    f->n = h->n;
    f->seg = malloc(sizeof(void*)*f->n);
    int i;
    for (i = 0; i < f->n; ++i) f->seg[i] = segment_alloc(h->seg[i]->n);
    return f;
}
void thread_hold_destroy(struct thread_hold *h)
{
    int i;
    for (i = 0; i < h->n; ++i) {
        struct segment *s = h->seg[i];
        free(s->counts);
        free(s);
    }
    free(h->seg);
    free(h);
}
struct BarcodeRegion {
    int rd; // read 1 or 2
    int start;
    int end;
    int dist;
    char **white_list;
    int len;
    int n_wl;
    strhash_t *wlhash;
    // It is NOT Thread safe. delete it at 11 April 2019
    // reset to 0 before counting. 
    // int q30_bases;
    // int bases;
    // int exact_match;
};

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
struct NameCountPair {
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

    // Cell barcodes found in white list, if no white list all barcodes treat as background
    struct NameCountPair *names;
    int n_name;
    int m_name;    
    strhash_t *cbhash;

    // Background
    struct NameCountPair *bgnames;
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
    FILE *report_fp; // old report handler
    FILE *html_report_fp;
    FILE *barcode_dis_fp;
    
    // hold thread safe data
    struct thread_hold **hold;

    // input file streaming
    struct fastq_handler *fastq;

    // stats of reads
    uint64_t raw_reads;
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
    .html_report_fp = NULL,
    .barcode_dis_fp = NULL,
    .hold = NULL,
    .fastq = NULL,

    .raw_reads = 0,
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
    const char *platform;
    const char *version;
    // consistant with white list if set
    const char *cell_barcode_tag; 
    const char *sample_barcode_tag;

    const char *raw_cell_barcode_tag;
    const char *raw_cell_barcode_qual_tag;
    int n_cell_barcode;
    struct BarcodeRegion *cell_barcodes;

    const char *raw_sample_barcode_tag;
    const char *raw_sample_barcode_qual_tag;
    int n_sample_barcode; // usually == 1, sometimes == 2,
    struct BarcodeRegion *sample_barcodes;

    const char *umi_tag;
    const char *umi_qual_tag;
    // UMI, usually random generated and located at one region
    struct BarcodeRegion *UMI;

    // clean read sequence
    struct BarcodeRegion *read_1; 
    struct BarcodeRegion *read_2;  // for single ends, read2 == NULL
} config = {
    .platform = NULL,
    .version = NULL,
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
           if (node->v.str) config.platform = strdup(node->v.str);
        }
        else if (strcmp(node->key, "version") == 0) {
            if (node->v.str) config.version = strdup(node->v.str);
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
            config.cell_barcodes = calloc(node->n,sizeof(struct BarcodeRegion));

            int j;
            for (j = 0; j < node->n; ++j) {
                const kson_node_t *n1 = kson_by_index(node, j);                
                if (n1 == NULL) error("cell barcode is empty.");
                if (n1->type != KSON_TYPE_BRACE) error("Format error. \"cell barcode\":[{},{}]");
                if (n1->n == 0) continue; // empty record
                struct BarcodeRegion *br = &config.cell_barcodes[j];
                
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
            config.read_1 = malloc(sizeof(struct BarcodeRegion));
            struct BarcodeRegion *br = config.read_1;
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
            config.read_2 = malloc(sizeof(struct BarcodeRegion));
            struct BarcodeRegion *br = config.read_2;
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
            config.UMI = malloc(sizeof(struct BarcodeRegion));
            struct BarcodeRegion *br = config.UMI;
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
        struct BarcodeRegion *br = &config.cell_barcodes[i];
        if (br->n_wl == 0) continue;
        br->wlhash = kh_init(str);
        int j;
        khint_t k;
        int ret;
        for (j = 0; j < br->n_wl; j++) {
            int len = strlen(br->white_list[j]);
            if (len != br->len) error("Inconsistance white list length. %d vs %d, %s", br->len, len, br->white_list[j]);
            if (br->dist>3) error("Set too much distance for cell barcode, allow 3 distance at max.");
            if (br->dist > len/2) error("Allow distance greater than half of barcode! Try to reduce distance.");
            k = kh_put(str, br->wlhash, br->white_list[j], &ret);
            if (!ret) error("Duplicated white list %s", br->white_list[j]);
            kh_val(br->wlhash, k) = j+1;
        }
    }
    kson_destroy(json);
}

struct seqlite {
    char *seq;
    char *qual;
};

void seqlite_destory(struct seqlite *s)
{
    free(s->seq);
    if (s->qual) free(s->qual);    
    free(s);
}

struct seqlite *extract_tag(struct bseq *b, const struct BarcodeRegion *r, struct BRstat *stat)
{
    if (b == NULL || r == NULL) return NULL;
    
    struct seqlite *p = malloc(sizeof(*p));

    // memset(stat, 0, sizeof(struct BRstat));
    
    char *s = NULL;
    char *q = NULL;
    if (r->rd == 1) {
        s = b->s0 + r->start -1;
        q = b->q0 ? b->q0 + r->start -1 : NULL;
    }
    else {
        s = b->s1 + r->start -1;
        q = b->q1 ? b->q1 + r->start -1 : NULL;
    }
    int l = r->end - r->start + 1;
    kstring_t seq = {0,0,0};
    kstring_t qual = {0,0,0};
    kputsn(s, l, &seq);
    if (q) kputsn(q, l, &qual);
    p->seq = seq.s;
    p->qual = qual.s;
    int i;
    for (i = 0; i < l; ++i) 
        if (p->qual[i]-33 >= 30) stat->q30_bases++;
    stat->bases = +l;
    return p;
}
// credit to https://github.com/wooorm/levenshtein.c
int levenshtein(char *a, char *b, int l) {
    int *cache = calloc(l, sizeof(int));
    int index = 0;
    int bIndex = 0;
    int distance;
    int bDistance;
    int result;
    char code;

    // initialize the vector.
    while (index < l) {
        cache[index] = index + 1;
        index++;
    }

    // Loop.
    while (bIndex < l) {
        code = b[bIndex];
        result = distance = bIndex++;
        index = 0;

        while (++index < l) {
            bDistance = code == a[index] ? distance : distance + 1;
            distance = cache[index];

            cache[index] = result = distance > result
                ? bDistance > result
                ? result + 1
                : bDistance
                : bDistance > distance
                ? distance + 1
                : bDistance;
        }
    }
    
    free(cache);    
    return result;
}

// -1 on unfound, 0 on No found on white list, >0 for iterater of white lists
int check_whitelist(char *s, const struct BarcodeRegion *r, int *exact_match)
{
    if (r->n_wl == 0) return 0;
    int len = strlen(s);
    if (len != r->len) error("Trying to check inconsistance length sequence.");
    
    khint_t k;
    k = kh_get(str, r->wlhash, s);
    if (k != kh_end(r->wlhash)) {
        *exact_match = 1;
        int id = kh_val(r->wlhash, k);
        return id;
    }
    if (r->dist > 0) {
        int ret = -1;
        int i;
        for (i = 0; i < r->n_wl; ++i) {
            char *f = r->white_list[i];
            int dist = levenshtein(f, s, len);
            if (dist <= r->dist) {
                if (ret == -1) ret = i+1;
                else return -1; // more than one barcode statisfied, skip it
            }
        }
        return ret;
    }
    return -1;
}
static void update_rname(struct bseq *b, const char *tag, char *s)
{
    kstring_t str = {0,0,0};
    kputs(b->n0, &str);
    kputs("|||", &str);
    kputs(tag, &str);
    kputs("|||", &str);
    kputs(s, &str);
    free(b->n0);
    b->n0 = str.s;
}

struct BRstat *extract_barcodes(struct bseq *b,
                                int n,
                                const struct BarcodeRegion *r,
                                const char *tag,
                                const char *raw_tag,
                                const char *raw_qual_tag,                                
                                const char *run_code,
                                struct thread_hold *hold
    )
{
    if (tag == NULL && raw_tag == NULL) return NULL;
    
    kstring_t str = {0,0,0};
    kstring_t qual = {0,0,0};
    kstring_t tag_str = {0,0,0};

    struct BRstat *stat = malloc(sizeof(struct BRstat));
    memset(stat, 0, sizeof(struct BRstat));
    stat->exact_match = 1;

    int i;
    for (i = 0; i < n; ++i) {

        const struct BarcodeRegion *br = &r[i];
        struct segment *seg = hold == NULL ? NULL : hold->seg[i];
        struct seqlite *s = extract_tag(b, br, stat);
        if (s == NULL) goto failed_check_barcode;

        int exact_match = 0;
        int ret = check_whitelist(s->seq, br, &exact_match);

        if (raw_tag) kputs(s->seq, &str);
        if (raw_qual_tag && s->qual) kputs(s->qual, &qual);
        
        if (tag) {
            if (ret>0)
                kputs(br->white_list[ret-1], &tag_str);
            else // in case no white list
                kputs(s->seq, &tag_str);
        }
        seqlite_destory(s);

        if (br->n_wl == 0)             
            continue;
        
        if (ret <= 0) {
            stat->filter = 1;
            continue;
        }
        if (ret > 0) {
            if (exact_match == 1) seg->counts[ret-1].matched++;
            else seg->counts[ret-1].corrected++;
        }
        
        if (ret == 0 || exact_match == 0) {
            stat->exact_match = 0;
            continue;
        }        
    }
    
    if (run_code) {
        kputc('-', &tag_str);
        kputs(run_code, &tag_str);
        struct fq_data *data = (struct fq_data*)b->data;
        data->bc_str = strdup(tag_str.s);
    }
    
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
                                 const struct BarcodeRegion *r,
                                 const char *tag,
                                 const char *raw_tag,
                                 const char *raw_qual_tag)
{
    return extract_barcodes(b, n, r, tag, raw_tag, raw_qual_tag, NULL, NULL);
}
struct BRstat *extract_cell_barcode_reads(struct bseq *b,
                               int n,
                               const struct BarcodeRegion *r,
                               const char *tag,
                               const char *raw_tag,
                               const char *raw_qual_tag,                               
                               const char *run_code,
                               struct thread_hold *hold)
{
    return extract_barcodes(b, n, r, tag, raw_tag, raw_qual_tag, run_code, hold);
}

struct BRstat *extract_umi(struct bseq *b, const struct BarcodeRegion *r, const char *tag, const char *qual_tag)
{
    if (tag == NULL) return NULL;
    struct BRstat *stat = malloc(sizeof(*stat));
    memset(stat, 0, sizeof(*stat));
    //int i;
    kstring_t str = {0,0,0};
    kstring_t qual = {0,0,0};
    struct seqlite *s = extract_tag(b, r, stat);
    if (qual_tag && s->qual)
        kputs(s->qual, &qual);
    if (tag)
        kputs(s->seq, &str);
    
    seqlite_destory(s);
   
    if (tag) update_rname(b, tag, str.s);
    if (qual_tag) update_rname(b, qual_tag, qual.s);

    if (str.m) free(str.s);
    if (qual.m) free(qual.s);
    return stat;
}
                
struct BRstat *extract_reads(struct bseq *b, const struct BarcodeRegion *r1, const struct BarcodeRegion *r2)
{
    assert(r1);
    struct BRstat *stat = malloc(sizeof(struct BRstat));
    memset(stat, 0, sizeof(struct BRstat));
    struct seqlite *s1 = extract_tag(b, r1, stat);
    struct seqlite *s2 = extract_tag(b, r2, stat);

    free(b->s0);
    if (b->q0) free(b->q0);
    b->s0 = strdup(s1->seq);
    if (s1->qual) b->q0 = strdup(s1->qual);

    if (b->s1) {
        free(b->s1);
        b->s1 = NULL;
        b->l1 = 0;
    }
    if (b->q1) {
        free(b->q1);
        b->q1 = NULL;
    }
    if (s2 && s2->seq) {
        b->s1 = strdup(s2->seq);
        b->l1 = strlen(s2->seq);
    }
    
    if (s2 && s2->qual) b->q1 = strdup(s2->qual);

    seqlite_destory(s1);
    if (s2) seqlite_destory(s2);
    
    return stat;
}

static void *run_it(void *_p, int idx)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    struct args *opts = p->opts;
    struct thread_hold *hold = opts->hold[idx];

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
                opts->run_code,
                hold
                );
            
            if (cell_stat == NULL || cell_stat->filter == 1) {
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

            if (opts->bgiseq_filter) {
                if (b->q0) {
                    int k;
                    int bad_bases = 0;
                    for (k = 0; k < 15 && k <b->l0; ++k) {
                        if (b->q0[k]-33<10) bad_bases++;
                    }
                    if (bad_bases >2) b->flag = FQ_FLAG_READ_QUAL;
                    else {
                        if (b->q1) {
                            for (k = 0; k < 15 && k <b->l1; ++k) {
                                if (b->q1[k]-33<10) bad_bases++;
                            }
                            if (bad_bases >2) b->flag = FQ_FLAG_READ_QUAL;
                        }
                    }                
                }
            }

            if (opts->qual_thres > 0 && b->q0) { 
                int k;
                int ave = 0;
                
                for (k = 0; k < b->l0; k++) ave+=b->q0[k]-33;
                if (ave/k < opts->qual_thres) {
                    b->flag = FQ_FLAG_READ_QUAL;
                    continue;
                }

                if (b->l1 > 0 && b->q1) {
                    ave = 0;
                    for (k = 0; k < b->l1; k++) ave+=b->q1[k]-33;
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
            fprintf(fp1, "%c%s\n%s\n", b->q0 ? '@' : '>', b->n0, b->s0);
            if (b->q0) fprintf(fp1, "+\n%s\n", b->q0);
            if (b->l1 > 0) {
                fprintf(fp2, "%c%s\n%s\n", b->q1 ? '@' : '>', b->n0, b->s1);
                if (b->q1) fprintf(fp2, "+\n%s\n", b->q1);
            }

            if (data->bc_str && opts->cbhash) {
                khint_t k;
                k = kh_get(str, opts->cbhash, (char*)data->bc_str);
                if (k == kh_end(opts->cbhash)) {
                    if (opts->n_name == opts->m_name) {
                        opts->m_name += 10000;
                        opts->names = realloc(opts->names,opts->m_name *sizeof(struct NameCountPair));                        
                    }
                    struct NameCountPair *pair = &opts->names[opts->n_name];
                    pair->name = strdup(data->bc_str);
                    pair->count = 1;
                    k = kh_put(str,opts->cbhash, pair->name, &ret);
                    kh_val(opts->cbhash, k) = opts->n_name;
                    opts->n_name++;
                }
                else {
                    int id = kh_val(opts->cbhash, k);
                    struct NameCountPair *pair = &opts->names[id];
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
                        opts->bgnames = realloc(opts->bgnames,opts->m_bg *sizeof(struct NameCountPair));                        
                    }
                    struct NameCountPair *pair = &opts->bgnames[opts->n_bg];
                    pair->name = strdup(data->bc_str);
                    pair->count = 1;
                    k = kh_put(str,opts->bghash, pair->name, &ret);
                    kh_val(opts->bghash, k) = opts->n_bg;
                    opts->n_bg++;
                }
                else {
                    int id = kh_val(opts->bghash, k);
                    struct NameCountPair *pair = &opts->bgnames[id];
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
        opts->barcode_exactly_matched += data->cr_exact_match;
        if (data->bc_str) free(data->bc_str);
        free(data);
    }
    bseq_pool_destroy(p);
}
int cmpfunc (const void *a, const void *b)
{    
    return ( ((struct NameCountPair*)b)->count - ((struct NameCountPair*)a)->count );
}
void cell_barcode_count_pair_write()
{
    if (args.n_name && args.cbhash) {
        qsort(args.names, args.n_name, sizeof(struct NameCountPair), cmpfunc);
        int i;
        for (i = 0; i < args.n_name; ++i) {
            fprintf(args.cbdis_fp, "%s\t%d\n", args.names[i].name, args.names[i].count);
            free(args.names[i].name);
        }
        free(args.names);
        fclose(args.cbdis_fp);
    }
    kh_destroy(str,args.cbhash);
}
void report_write()
{
    fprintf(args.report_fp, "[\n");
    fprintf(args.report_fp, "\t{\"name\":\"Number of Fragments\", \"value\":\"%"PRIu64"\"},\n", args.raw_reads);
    fprintf(args.report_fp, "\t{\"name\":\"Fragments with Exactly Matched Barcodes\", \"value\":\"%"PRIu64"\"},\n", args.barcode_exactly_matched);
    fprintf(args.report_fp, "\t{\"name\":\"Fragments with Failed Barcodes\", \"value\":\"%"PRIu64"\"},\n", args.filtered_by_barcode);
    fprintf(args.report_fp, "\t{\"name\":\"Fragments Filtered on Low Qulity\", \"value\":\"%"PRIu64"\"},\n", args.filtered_by_lowqual);
    fprintf(args.report_fp, "\t{\"name\":\"Fragments Filtered on Unknown Sample Barcodes\", \"value\":\"%"PRIu64"\"},\n", args.filtered_by_sample);
    fprintf(args.report_fp, "\t{\"name\":\"Q30 bases in Cell Barcode\", \"value\":\"%.1f%%\"},\n",args.bases_cell_barcode == 0 ? 0 : (float)args.q30_bases_cell_barcode/(args.bases_cell_barcode+1)*100);
    fprintf(args.report_fp, "\t{\"name\":\"Q30 bases in Sample Barcode\", \"value\":\"%.1f%%\"},\n", args.bases_sample_barcode == 0 ? 0 : (float)args.q30_bases_sample_barcode/(args.bases_sample_barcode+1)*100);
    fprintf(args.report_fp, "\t{\"name\":\"Q30 bases in UMI\", \"value\":\"%.1f%%\"},\n", args.bases_umi == 0 ? 0 : (float)args.q30_bases_umi/(args.bases_umi+1)*100);
    fprintf(args.report_fp, "\t{\"name\":\"Q30 bases in Reads\", \"value\":\"%.1f%%\"}\n", (float)args.q30_bases_reads/(args.bases_reads+1)*100);
    fprintf(args.report_fp, "]\n");
    if (args.report_fname) fclose(args.report_fp);
}
void full_details()
{
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
        struct BarcodeRegion *br = &config.cell_barcodes[i];
        struct segment *s0 = args.hold[0]->seg[i];
        fprintf(args.barcode_dis_fp, "# Read %d, %d-%d\n", br->rd, br->start, br->end);
        for (j = 0; j < br->n_wl; ++j)
            fprintf(args.barcode_dis_fp, "%s\t%"PRIu64"\t%"PRIu64"\n", br->white_list[j], s0->counts[j].matched, s0->counts[j].corrected);
    }    
}
static void memory_release()
{
    if (args.r1_fp) gzclose(args.r1_fp);
    if (args.r2_fp) gzclose(args.r2_fp);
    if (args.out1_fp) fclose(args.out1_fp);
    if (args.out2_fp) fclose(args.out2_fp);
    fclose(args.barcode_dis_fp);
    fastq_handler_destory(args.fastq);
    int i;
    for (i = 0; i < args.n_thread; ++i) thread_hold_destroy(args.hold[i]);
    free(args.hold);    
}
static int parse_args(int argc, char **argv)
{
    if ( argc == 1 ) return usage();

    int i;
    const char *thread = NULL;
    const char *chunk_size = NULL;    
    const char *qual_thres = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return usage();

        if (strcmp(a, "-1") == 0) var = &args.out1_fname;
        else if (strcmp(a, "-2") == 0) var = &args.out2_fname;
        else if (strcmp(a, "-config") == 0) var = &args.config_fname;
        else if (strcmp(a, "-cbdis") == 0) var = &args.cbdis_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-r") == 0) var = &chunk_size;
        else if (strcmp(a, "-run") == 0) var = &args.run_code;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-dis") == 0) var = &args.dis_fname;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;
        else if (strcmp(a, "-f") == 0) {
            args.bgiseq_filter = 1;
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
    if (args.cbdis_fname == NULL) error("-cbdis is required.");
    config_init(args.config_fname);
        
    if (thread) args.n_thread = str2int((char*)thread);
    if (chunk_size) args.chunk_size = str2int((char*)chunk_size);
    assert(args.n_thread >= 1 && args.chunk_size >= 1);
    if (qual_thres) {
        args.qual_thres = str2int((char*)qual_thres);
        LOG_print("Average quality below %d will be drop.", args.qual_thres);
    }


    args.hold = malloc(sizeof(void*)*args.n_thread);
    memset(args.hold, 0, sizeof(void*)*args.n_thread);
    
    if (config.n_cell_barcode > 0) { // segments        
        struct thread_hold *h1 = malloc(sizeof(struct thread_hold));
        h1->n = config.n_cell_barcode;
        h1->seg = malloc(sizeof(void*)*h1->n);
        
        for (i = 0; i < config.n_cell_barcode; ++i) { // barcodes per segment
            int wl = config.cell_barcodes[i].n_wl;
            h1->seg[i] = malloc(sizeof(struct segment));
            h1->seg[i]->n = wl;
            h1->seg[i]->counts = malloc(wl*sizeof(struct bcount));
            memset(h1->seg[i]->counts, 0, sizeof(struct bcount)*wl);
        }

        args.hold[0] = h1;
        for (i = 1; i < args.n_thread; ++i) args.hold[i] = thread_hold_fork(h1);
    }

    if (args.r1_fname == NULL && (!isatty(fileno(stdin)))) args.r1_fname = "-";
    if (args.r1_fname == NULL) error("Fastq file(s) must be set.");
        
    if (args.report_fname) {
        args.report_fp = fopen(args.report_fname, "w");
        CHECK_EMPTY(args.report_fp, "%s : %s.", args.report_fname, strerror(errno));
    }
    else error("-report must be set.");

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
    }   
        
    if (args.run_code == NULL) args.run_code = strdup("1");

    if (args.dis_fname != NULL) {
        args.barcode_dis_fp = fopen(args.dis_fname, "w");
        CHECK_EMPTY(args.barcode_dis_fp, "%s : %s.", args.dis_fname, strerror(errno));
    }
    else {
        args.barcode_dis_fp = stderr;
    }
    
    args.cbhash = kh_init(str);
    args.bghash = kh_init(str);
    
    return 0;
}

int fastq_prase_barcodes(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return 1;

    int nt = args.n_thread;
    
    struct thread_pool *p = thread_pool_init(nt);
    struct thread_pool_process *q = thread_pool_process_init(p, nt*2, 0);
    struct thread_pool_result *r;

    for (;;) {
        struct bseq_pool *b = bseq_read(args.fastq, &args);
        if (b == NULL) break;
        int block;
        do {
            block = thread_pool_dispatch2(p, q, run_it, b, 0);
            if ((r = thread_pool_next_result(q))) {
                struct bseq_pool *d = (struct bseq_pool *)r->data;
                write_out(d);
            }
            thread_pool_delete_result(r, 0);
        }
        while (block == -1);
    }
    thread_pool_process_flush(q);

    while ((r = thread_pool_next_result(q))) {
        struct bseq_pool *d = (struct bseq_pool*)r->data;
        write_out(d);
        thread_pool_delete_result(r, 0);
    }

    thread_pool_process_destroy(q);
    thread_pool_destroy(p);

    cell_barcode_count_pair_write();

    report_write();
    // todo: html report
    full_details();
    
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    
    return 0;
}
