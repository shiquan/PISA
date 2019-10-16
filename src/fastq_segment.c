/*
  Parse and trim single cell LFR adaptors and barcodes
  
  demo:
{
    "pattern":"CCCCAGAACGACATGGCTACGAAGTCGGANNNNNNNNNNNNNNNNNNNCCTTATCAGCGTCCGACTTCGTAGCCATGTCGTTCTGCG(T1)CCTTCC(T2)CGATG(UM)[polyT]",
    "segments":[
        { // LFR cell barcode
            "tag":"CB", 
            "seg":[  // two segments
                {
                    "tag":"T1",
                    "length":"10",
                },
                {
                    "tag":"T2",
                    "length":"10",
                }
            ]
        },
        { // UMI
            "tag":"UM",
            "length":"10",
        }
    ]
}

 */
#include "utils.h"
#include "json_config.h"
#include "dict.h"
#include "kson.h"
#include "number.h"
#include "dict.h"
#include "read_thread.h"
#include "read_tags.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/thread_pool.h"

struct ref;

static struct args {
    const char *input_fname;
    const char *config_fname;
    const char *output_fname;
    int seed_length;    
    struct ref *r;
    
    const char *phase_tag; // reads have same RB tag in each block from one template
    int n_thread;
    int input_pe;
    int filter_reads;

    int keep_all;
    struct dict *tag_dict;
    FILE *out;    
} args = {
    .input_fname = NULL,
    .config_fname = NULL,
    .output_fname = NULL,
    .seed_length = 6,
    
    .r = NULL,
    .keep_all = 0,
    .tag_dict = NULL,
    .phase_tag = NULL,    
    .filter_reads = 0,
    .n_thread = 1,
    .input_pe = 0,
    .out = NULL,
};

extern char* check_circle(char *seq);

#define MIN_SEED_LENGTH    6
#define MAX_SEED_LENGTH    20
struct hit {
    int loc;
    int strand; // 0 for plus, -1 for minus
};
struct hits {
    int n;
    struct hit *hit;
};

KHASH_MAP_INIT_STR(hit, struct hits*)
KHASH_MAP_INIT_STR(str, int)

#define MAX_MISMATCH 3

struct segment1 {
    int idx;
    char *tag;
    int l;
    int dist;
    int n;
    char **wl; // white list
    kh_str_t *hash; // white list hash;    
};

struct segment {
    char *tag;
    int n;
    struct segment1 *s;
};

struct tag {
    int i1;
    int i2;
};

KHASH_MAP_INIT_STR(tag, struct tag)

struct ref_pat {

    struct tag **tags; // same length with seq
    char *seq; // ACGTN, N for tag, B for polyA, P for polyT
    int l_seq;

    int n_tag; // all tags for alloc
    
    int n;
    struct segment *segs;
};
struct ref {
    struct ref_pat *ref;
    struct ref_pat *rev;
    int n, m;
    char **kmers;
    kh_hit_t *map;
    kh_tag_t *hash;
};

static struct ref_pat *ref_pat_alloc()
{
    struct ref_pat *r = malloc(sizeof(*r));
    memset(r, 0, sizeof(*r));
    return r;
}
static void ref_pat_destroy(struct ref_pat *r)
{
    int i;
    for (i = 0; i < r->n; ++i) {
        struct segment *s = &r->segs[i];
        if (s->tag) free(s->tag);
        int k;
        for (k = 0; k < s->n; ++k) {
            struct segment1 *s0 = &s->s[k];
            if (s0->tag) free(s0->tag);
            if (s0->n > 0) {
                int j;
                for (j = 0; j < s0->n; ++j) free(s0->wl[j]);
                free(s0->wl);
                kh_destroy(str, s0->hash);
            }
        }
    }
    free(r->segs);
    free(r->seq);
    free(r->tags);
    free(r);
}

static void ref_destroy(struct ref *r)
{
    ref_pat_destroy(r->ref);
    ref_pat_destroy(r->rev);
    int i;
    for (i = 0; i < r->n; ++i) free(r->kmers[i]);
    free(r->kmers);
    kh_destroy(hit, r->map);
    kh_destroy(tag, r->hash);
    free(r);
}
static char *rev_seq(char *s, int l)
{
    char *r = malloc((l+1)*sizeof(char));
    int i;
    for (i = 0; i < l; ++i) {
        switch(s[i]) {
            case 'A':
                r[l-i-1] = 'T';
                break;
            case 'C':
                r[l-i-1] = 'G';
                break;
            case 'G':
                r[l-i-1] = 'C';
                break;
            case 'T':
                r[l-i-1] = 'A';
                break;
            case 'N':
                r[l-i-1] = 'N';
                break;
            case 'B':
                r[l-i-1] = 'P';
                break;
            case 'P':
                r[l-i-1] = 'B';
                break;
                
            default:
                error("Try to reverse a non DNA sequence? \"%s\"", s);
        }
    }
    r[l] = '\0';
    return r;
}
static struct ref_pat *build_rev_pat(struct ref_pat *r)
{
    struct ref_pat *v = ref_pat_alloc();
    int i;
    v->l_seq = r->l_seq;
    v->tags = malloc(sizeof(void*)*r->l_seq);
    for (i = 0; i < r->l_seq; ++i) v->tags[i] = r->tags[r->l_seq-i-1];
    v->seq = rev_seq(r->seq, r->l_seq);
    v->n = r->n;
    v->n_tag = r->n_tag;
    v->segs = malloc(v->n*sizeof(struct segment));
    for (i = 0; i < v->n; ++i) {
        struct segment *ori = &r->segs[i];
        struct segment *copy = &v->segs[i];
        if (ori->tag) copy->tag = strdup(ori->tag);
        copy->n = ori->n;
        copy->s = malloc(ori->n*sizeof(struct segment1));
        int k;
        for (k = 0; k < ori->n; ++k) {
            struct segment1 *s1 = &ori->s[k];
            struct segment1 *s2 = &copy->s[k];
            memset(s2, 0, sizeof(*s2));
            s2->l = s1->l;
            s2->dist = s1->dist;
            s2->idx = s1->idx;
            if (s1->tag) s2->tag = strdup(s1->tag);
            if (s1->n) {
                s2->n = s1->n;
                s2->hash = kh_init(str);
                s2->wl = malloc(sizeof(char*)*s2->n);
                int j;
                khint_t ik;
                int ret;
                for (j = 0; j < s1->n; ++j) {
                    s2->wl[j] = rev_seq(s1->wl[j], s1->l);
                    ik = kh_put(str, s2->hash, s2->wl[j], &ret);
                    kh_val(s2->hash, ik) = j;
                }
            }
        }
    }
    return v;
}

static void build_kmers_ref(struct ref *r)
{
    char *ref[2];
    ref[0] = r->ref->seq;
    ref[1] = r->rev->seq;
    
    r->map = kh_init(hit);
    
    int strand = 0;
    for (; strand < 2; ++strand) {
        int i;
        int l = r->ref->l_seq;
        // debug_print("ref : %d %s", strand, ref[strand]);
        for (i = 0; i < l; ) {
            char *s = ref[strand];
            int j;
            kstring_t str = {0,0,0};
            if (s[i] == 'N') { i++; continue; }
            else if (s[i] == 'B')
                for (j = 0; j < args.seed_length; ++j) kputc('A', &str);
            else if (s[i] == 'P')
                for (j = 0; j < args.seed_length; ++j) kputc('T', &str);
            else {
                for (j = 0; j < args.seed_length && i+j < r->ref->l_seq; j++)
                    if (s[i+j] == 'N' || s[i+j] == 'B' || s[i+j] == 'P') break;
                if (j < args.seed_length) { i+=j; continue; }
                kputsn(s+i, args.seed_length, &str);
            }
            kputs("", &str);
            khint_t k;
            int ret;
            k = kh_put(hit, r->map, str.s, &ret);
            if (ret) {
                struct hits *h = malloc(sizeof(struct hits));
                h->n = 1;
                h->hit = malloc(sizeof(struct hit));
                h->hit->strand = strand;
                h->hit->loc = i;
                kh_val(r->map, k) = h;
                if (r->n == r->m) {
                    r->m += 100;
                    r->kmers = realloc(r->kmers, r->m*sizeof(char*));
                }
                r->kmers[r->n++] = str.s; // keep key
                // debug_print("kmers : %s, %d, %d", str.s, i, strand);
            }
            else {
                struct hits *h = kh_val(r->map, k);
                h->hit = realloc(h->hit, (h->n+1)*sizeof(struct hit));
                struct hit *hh = &h->hit[h->n++];
                hh->strand = strand;
                hh->loc = i;
                free(str.s);
            }            
            i++;
        }
    }
}

static int parse_seg1(struct segment1 *s, const kson_node_t *n)
{
    if (n == NULL) error("Segment record is empty.");
    int i;
    for (i = 0; i < n->n; ++i) {
        const kson_node_t *n1 = kson_by_index(n, i);
        if (n1 == NULL) continue;
        if (n1->key == NULL) error("Format error.");
        if (strcmp(n1->key,"tag") == 0) {            
            s->tag = strdup(n1->v.str);
        }
        else if (strcmp(n1->key, "length") == 0) {
            s->l = str2int(n1->v.str);
        }
        else if (strcmp(n1->key, "distance") == 0) {
            s->dist = str2int(n1->v.str);
        }
        else if (strcmp(n1->key, "white list") == 0) {
            if (n1->type != KSON_TYPE_BRACKET) error("Format error.\"white list\":[]");
            s->n = n1->n;
            s->wl = malloc(s->n*sizeof(char*));
            s->hash = kh_init(str);
            int m;
            for (m = 0; m <n1->n; ++m) {
                const kson_node_t *n2 = kson_by_index(n1,m);
                if (n2->v.str == NULL) error ("Empty record at white list.");
                s->wl[m] = strdup(n2->v.str);
                khint_t k;
                int ret;
                k = kh_put(str, s->hash, s->wl[m], &ret);
                if (!ret) error ("Duplicate record at white list. %s", s->wl[m]);
                kh_val(s->hash, k) = m;
            }
        }            
        else error("Unknown key %s", n->key);
    }
    
    return 0;
}
static int parse_seg(struct segment *s, const kson_node_t *n0)
{
    if (n0 == NULL) error("Segment record is empty.");
    memset(s, 0, sizeof(*s));
    s->n = 1;
    s->s = malloc(sizeof(struct segment1));
    struct segment1 *s0 = s->s;
    memset(s0, 0, sizeof(struct segment1));
    s0->tag = NULL;
    int is_collect = 0;
    int k;
    for (k = 0; k < n0->n; ++k) {
        const kson_node_t *n = kson_by_index(n0, k);
        if (strcmp(n->key, "tag") == 0) {
            s->tag = strdup(n->v.str);
        }
        else if (strcmp(n->key, "seg") == 0) {
            if (n->type != KSON_TYPE_BRACKET) error("Format error. \"seg\":[]");
            s->n = n->n;
            s->s = realloc(s->s, s->n*sizeof(struct segment1));
            int i;
            for (i = 0; i < s->n; ++i) {
                const kson_node_t *n1 = kson_by_index(n, i);            
                memset(&s->s[i], 0, sizeof(struct segment1));
                parse_seg1(&s->s[i], n1);
            }
            is_collect = 1;
        }
        else if (strcmp(n->key, "length") == 0) {
            if (is_collect == 1) error("Unknown format. \"seg\" is conflict with \"length\".");
            s0->l = str2int(n->v.str);
        }
        else if (strcmp(n->key, "distance") == 0) {
            if (is_collect == 1) error("Unknown format. \"seg\" is conflict with \"distance\".");
            s0->dist = str2int(n->v.str);
        }
        else if (strcmp(n->key, "white list") == 0) {
            if (is_collect == 1) error("Unknown format. \"seg\" is conflict with \"white list\".");
            if (n->type != KSON_TYPE_BRACKET) error("Format error.\"white list\":[]");
            s0->n = n->n;
            s0->wl = malloc(s0->n*sizeof(char*));
            s0->hash = kh_init(str);
            int m;
            for (m = 0; m <n->n; ++m) {
                const kson_node_t *n2 = kson_by_index(n,m);
                if (n2->v.str == NULL) error ("Empty record at white list.");
                s0->wl[m] = strdup(n2->v.str);
                khint_t k;
                int ret;
                k = kh_put(str, s0->hash, s0->wl[m], &ret);
                if (!ret) error ("Duplicate record at white list. %s", s0->wl[m]);
                kh_val(s0->hash, k) = m;                                
            }        
        }
        else error("Unknown key %s", n->key);
    }
    return 0;
}
static struct ref_pat *config_parse(const char *fn)
{
    char *config_str = json_config_open(fn);
    if (config_str == NULL) error("Empty configure file.");
    kson_t *json = kson_parse(config_str);
    free(config_str);
    
    struct ref_pat *ref = ref_pat_alloc();

    const kson_node_t *root = json->root;
    if (root == NULL) error("Format error. Root node is emtpy.");
    int i;
    for (i = 0; i < root->n; ++i) {
        const kson_node_t *node = kson_by_index(root, i);
        if (node == NULL) continue;
        if (node->key == NULL) error("Format error. Node key is empty.");
        if (strcmp(node->key, "pattern") == 0) {
            ref->seq = strdup(node->v.str);
        }
        else if (strcmp(node->key, "segments") == 0) {

            if (node->type != KSON_TYPE_BRACKET) error("Format error.\"segments\":[{},{}]");

            ref->n = node->n;
            if (ref->n == 0) error("No segment records in the config file.");
            ref->segs = malloc(ref->n*sizeof(struct segment));
            int j;
            for (j =0; j <node->n; ++j) {
                const kson_node_t *n = kson_by_index(node, j);
                if (n == NULL) error("Segment record is empty.");
                if (n->type!= KSON_TYPE_BRACE) error("Format error.\"segments\":[{},{}]");
                struct segment *seg = &ref->segs[j];
                memset(seg, 0, sizeof(struct segment));
                parse_seg(seg, n);
            }
        }
        else error("Unknown key, %s", node->key);
    }
    kson_destroy(json);
    
    if (ref->seq == NULL) error("\"pattern\" must be set in the configure file.");
    return ref;
}

static void config_update_tag_index(struct ref *r)
{
    int i;
    int n;
    struct ref_pat *ref = r->ref;
    for (i = 0, n = 0; i < ref->n; ++i) 
        n += ref->segs[i].n;

    if (n == 0) error("Empty config.");
    r->hash = kh_init(tag);
    int idx = 0;
    for (i = 0; i < ref->n; ++i) {
        khint_t k;
        int ret;
        struct segment *s = &ref->segs[i];
        if (s->tag == NULL) error("No tag found at segment.");
        if (s->n == 1) {
            k = kh_put(tag, r->hash, s->tag, &ret);
            if (!ret) error("Duplicate tag %s", s->tag);
            struct tag *t = &kh_val(r->hash, k);
            t->i1 = i; t->i2 = 0;
            s->s[0].idx = idx++;
        }
        else {
            int j;
            for (j = 0; j < s->n; ++j) {
                struct segment1 *s1 = &s->s[j];
                if (s1->tag == NULL) error("No tag found at segment.");
                k = kh_put(tag, r->hash, s1->tag, &ret);
                if (!ret) error("Duplicate tag %s", s1->tag);
                struct tag *t = &kh_val(r->hash, k);
                t->i1 = i; t->i2 = j;
                s1->idx = idx++;
            }
        }
    }
    ref->n_tag = idx;
}

static void build_ref_pattern(struct ref *ref)
{
    struct ref_pat *r = ref->ref;
    char *pat = r->seq;
    kstring_t str = {0,0,0};
    kstring_t seq = {0,0,0};
    int is_tag = 0;
    int i;
    int l = strlen(pat);
    r->tags = malloc(l*sizeof(void*));
    memset(r->tags, 0, sizeof(void*)*l);
        
    for (i = 0; i < l; ++i) {
        switch(pat[i]) {
            case '\0':
                do {
                    if (is_tag == 1) error("Unkown format.");
                    kputsn(str.s, str.l, &seq);
                } while(0);
                break;
                
            case '(':
                do {
                    if (str.l != 0) kputsn(str.s, str.l, &seq);
                    str.l = 0;
                    if (is_tag == 1) error("Unknown format.");
                    is_tag = 1;
                } while(0);
                break;
                    
            case ')':
                do {
                    if (is_tag == 0) error("Unknown format.");
                    if (str.l != 2) error("Tag name is too long. Only support 2 character, such 'AS','CB' etc.");
                        
                    khint_t k;
                    k = kh_get(tag, ref->hash, str.s);
                    if (k == kh_end(ref->hash)) error("No tag %s found at configure.", str.s);
                    struct tag *t = &kh_val(ref->hash, k);
                    r->tags[seq.l] = t;
                    kputc('N', &seq);
                    str.l = 0;
                    is_tag = 0;
                } while (0);
                break;

            case '[':
                do {
                    if (str.l != 0) kputsn(str.s, str.l, &seq);
                    str.l = 0;
                    if (is_tag == 1) error("Unknown format.");
                    is_tag = 1;
                } while(0);
                break;

            case ']':
                do {
                    if (is_tag == 0 || str.l == 0) error("Unknown format.");
                    if (strcmp(str.s, "polyA") == 0) kputc('B', &seq);
                    else if (strcmp(str.s, "polyT") == 0) kputc('P', &seq);
                    else error("Unknow tag, only support [polyA] or [polyT] for now.");
                    is_tag = 0;
                    str.l = 0;
                } while(0);
                break;
                
            case 'a': case 'A':
            case 'c': case 'C':                
            case 'g': case 'G':
            case 't': case 'T':
            case 'N':
                kputc(pat[i], &str);
                break;
                        
            default:
                do {
                    if (is_tag == 0) error("Unsupport DNA sequence in the pattern string. '%c'", pat[i]);
                    kputc(pat[i], &str);
                    // if (str.l > 2) error("Tag name is too long. Only support 2 character, such 'AS','CB' etc.");
                } while(0);
                break;

        }
    }

    if (is_tag == 1) error("Unknown format.");
    if (str.l) kputs(str.s, &seq);
    free(str.s);

    r->seq = seq.s;
    r->l_seq = seq.l;

    for (i = 0; i < r->l_seq; ++i) {
        switch(r->seq[i]) {
            case 'a':
                r->seq[i] = 'A';
                break;

            case 'c':
                r->seq[i] = 'C';
                break;
                
            case 'g':
                r->seq[i] = 'G';
                break;

            case 't':
                r->seq[i] = 'T';
                break;

            default:
                break;
        }
    }
    free(pat);
}

static struct ref *config_init(const char *fn)
{
    struct ref *r = malloc(sizeof(*r));
    memset(r, 0, sizeof(*r));
    r->ref = config_parse(fn);
    
    config_update_tag_index(r);
    
    build_ref_pattern(r);
    
    r->rev = build_rev_pat(r->ref);

    build_kmers_ref(r);

    return r;
}

static int usage()
{
    fprintf(stderr, "* Pick pre-designed segments from untigs.\n");
    fprintf(stderr, "Segment [options] in.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "-config [json]          Configure file.\n");    
    fprintf(stderr, "-o      [fastq]         Trimed fastq.\n");
    fprintf(stderr, "-sl     [INT]           Seed length for mapping consensus sequence.\n");
    fprintf(stderr, "-t      [INT]           Threads.\n");
    fprintf(stderr, "-tags   [CB,UMI]        Tags for each read block.\n");
    fprintf(stderr, "-pb     [PB]            Phase block tag.\n");
    fprintf(stderr, "-k                      Keep all reads even no segments detected.\n");
    fprintf(stderr, "-sum    summary.csv     Summary report.\n");
    return 1;
}

struct tag_val {
    int n;
    struct dict *tag;
    struct dict **dict;
};

static int check_pattern_right(char *s, const int start,  struct ref_pat *r, const struct hit *h, char **pat, int is_circle, int *checked)
{
    int len = strlen(s);
    int st1, st2;
    int mis = 0;
    int n_base = 0;
    *checked = 0;
    st1 = h->loc;
    if (r->seq[st1] != 'B' && r->seq[st1] != 'P') st1 += args.seed_length; // skip seed if not polyT/As
    st2 = start + args.seed_length;
    for ( ; st1 < r->l_seq;) {
        if (st2 >= len) { // adjust position for circle sequence
            if (is_circle) st2 = st2-len; // circle
            else { st2 = len; break; }
        }

        if (r->seq[st1] == 'N') {
            if (r->tags[st1] == NULL) {
                st1++; st2++;
            }
            else {
                struct tag *t = r->tags[st1];
                struct segment1 *g = &r->segs[t->i1].s[t->i2];
                kstring_t str = {0,0,0};
                if (len - st2 < g->l) {
                    if (is_circle) {
                        kputs(s+st2, &str);
                        kputsn(s, g->l+st2-len, &str);
                        kputs("", &str);
                    }
                    else return len; // return the last position
                }
                else {
                    kputsn(s+st2, g->l, &str); kputs("", &str);
                }
                if (g->n) {
                    khint_t k;
                    k = kh_get(str, g->hash, str.s);
                    if (k == kh_end(g->hash)) {
                        // todo: correct
                    }                    
                }
                if (pat[g->idx] != NULL) {
                    free(pat[g->idx]); // for polyA/T ends, may introduce bias if last base at tag is A or T
                    warnings("Double tag found.");
                }
                pat[g->idx] = str.s; 
                st2 += g->l;
                st1++;
            }
        }
        else if (r->seq[st1] == 'B') {
            // if (s[st2] == 'A') { st2++; continue;}
            //else st1++;
            st1++;
        }
        else if (r->seq[st1] == 'P') {
            if (s[st2] == 'T') { st2++; continue;}
            else  st1++;
        }
        else {
            if (r->seq[st1] != s[st2]) mis++;
            if (mis > MAX_MISMATCH) return -1;
            st1++, st2++;
            n_base++;
        }
    }
    // debug_print("right: %d,%d", n_base, mis);
    // todo: improve filtering
    if (n_base <= 30 && mis >= MAX_MISMATCH) return -1;
    if (n_base > 30 && (float)mis/n_base > 0.1) return -1;
    if (n_base > 5) *checked = 1; // pattern checked, not just seed
    return st2;
}
static int check_pattern_left(char *s, const int start, struct ref_pat *r, const struct hit *h, char **pat, int is_circle, int *checked)
{
    int st1, st2;
    int len = strlen(s);
    int mis = 0;
    int n_base = 0;
    *checked = 0;
    for (st1 = h->loc, st2= start; st1 >= 0;) {
        if (st2 < 0) {
            if (is_circle) st2 = len+st2;
            else break;
        }
        
        if (r->seq[st1] == 'N') {
            if (r->tags[st1] == NULL) {
                st1--; st2--;
            }
            else {
                struct tag *t = r->tags[st1];
                struct segment1 *g = &r->segs[t->i1].s[t->i2];
                kstring_t str = {0,0,0};
                if (st2 < g->l) {
                    if (is_circle) {
                        char *ss = s+len+st2-g->l;
                        kputsn(ss, g->l-st2, &str);
                        kputsn(s, st2, &str);
                        kputs("",&str);
                    }
                    else return 0; // truncated tag, ignore this tag.
                }
                else {
                    char *ss = s + st2 - g->l + 1;
                    kputsn(ss, g->l, &str);
                    kputs("", &str);
                }
                
                if (g->n) {
                    khint_t k;
                    k = kh_get(str, g->hash, str.s);
                    if (k == kh_end(g->hash)) {
                        // todo: correct
                    }                    
                }
                if (pat[g->idx] != NULL) {
                    free(pat[g->idx]);
                    warnings("Double tag found.");
                }
                pat[g->idx] = str.s;
                st2 -= g->l;
                st1--;
            }
        }
        else if (r->seq[st1] == 'B') { // TODO: for polyT/As. should redesigned here.
            // todo: improve performance here
            if (s[st2] == 'A') { st2--; continue;}
            else st1--;
        }
        else if (r->seq[st1] == 'P') {
            st1--;
        }
        else {
            if (r->seq[st1] != s[st2]) mis++;
            if (mis > MAX_MISMATCH) return -1;
            st1--, st2--;
            n_base++;
        }
    }    // todo: improve filtering
    if (n_base <= 30 && mis >= MAX_MISMATCH) return -1;
    if (n_base > 30 && (float)mis/n_base > 0.1) return -1;
    if (n_base > 5) *checked = 1;
    return st2+1;
}

static void pat2tagval(struct ref_pat *r, char **pat, int strand, struct tag_val *v)
{
    int i;
    for (i = 0; i < r->n; ++i) {
        struct segment *seg = &r->segs[i];
        if (seg->n > 1) {
            int j;        
            for (j = 0; j < seg->n; ++j) {
                char *a = pat[seg->s[j].idx];
                if (a == NULL) continue;
                int idx = dict_query(v->tag, seg->s[j].tag); // map tag to id
                assert(idx>=0);
                if (v->dict[idx] == NULL) v->dict[idx] = dict_init();
                if (strand == 1) {
                    char *rev = rev_seq(a, strlen(a));
                    dict_push(v->dict[idx], rev);
                    free(rev);
                }
                else {
                    dict_push(v->dict[idx], a);
                }
            }
        }
        else {
            char *a = pat[seg->s[0].idx];
            if (a == NULL) continue;
            int idx = dict_query(v->tag, seg->tag); // map tag to id
            assert(idx>=0);
            if (v->dict[idx] == NULL) v->dict[idx] = dict_init();
            if (strand == 1) {
                char *rev = rev_seq(a, strlen(a));
                dict_push(v->dict[idx], rev);
                free(rev);
            }
            else {
                dict_push(v->dict[idx], a);
            }
        }
    }
}
static char * tagval2name(struct ref *ref, struct tag_val *v)
{
    struct ref_pat *r = ref->ref;
    kstring_t str={0,0,0};
    int i;
    for (i = 0; i < r->n; ++i) {
        struct segment *seg = &r->segs[i];
        ksprintf(&str, "|||%s:Z:", seg->tag);
        if (seg->n == 1) {
            int idx = dict_query(v->tag, seg->tag);
            assert(idx>=0);
            if (v->dict[idx] == NULL) goto empty_bc;
            char *a;
            if (dict_size(v->dict[idx]) > 1) 
                a =dict_most_likely_key(v->dict[idx]);
            else
                a = dict_name(v->dict[idx], 0);
            if (a == NULL) goto empty_bc;
            kputs(a, &str);
        }
        else {
            int k;
            for (k = 0; k < seg->n; ++k) {
                int idx = dict_query(v->tag, seg->s[k].tag);            
                assert(idx>=0);
                if (v->dict[idx] == NULL) goto empty_bc;
                char *a;
                if (dict_size(v->dict[idx]) > 1)
                    a =dict_most_likely_key(v->dict[idx]);
                else
                    a = dict_name(v->dict[idx], 0);
                if (a==NULL) goto empty_bc;
                kputs(a, &str);
            }
        }
    }

    return str.s;
    
  empty_bc:
    if (str.m != 0) free(str.s);
    return NULL;
}

static char **check_pattern(char *s, int start, struct ref *ref, struct hit *h, int is_circle, struct tag_val *v, int *n_reads)
{
    struct ref_pat *r = h->strand == 0 ? ref->ref : ref->rev;
    *n_reads = 0;
    char **pat = malloc(r->n_tag*sizeof(void*));
    memset(pat, 0, r->n_tag*sizeof(void*));
    int l;
    int i;
    int s1, s2;
    l = strlen(s);
    int c1, c2;
    s1 = check_pattern_left(s, start, r, h, pat, is_circle, &c1);
    s2 = check_pattern_right(s, start, r, h, pat, is_circle, &c2);
    
    if (s1 == -1 || s2 == -1 || (c1==0 && c2==0)) {
        for (i = 0; i < r->n_tag; ++i)
            if (pat[i] != NULL) free(pat[i]);
        free(pat);
        return NULL;
    }

    pat2tagval(r, pat, h->strand, v);

    for (i = 0; i < r->n_tag; ++i)
        if (pat[i]) free(pat[i]);
    free(pat);

    *n_reads = 0;

    int n = 0;

    if (s1 == 0) {
        if (s2 == l) return NULL;
        *n_reads = 1;
        kstring_t str = {0,0,0};
        kputs(s+s2, &str);
        *n_reads = 1;
        char **ret = malloc(sizeof(void*)*(*n_reads));
        ret[0] = str.s;
        return ret;
    }
    else if (s2 == l) {
        kstring_t str = {0,0,0};
        kputsn(s, s1, &str);
        kputs("", &str);
        *n_reads = 1;
        char **ret = malloc(sizeof(void*)*(*n_reads));
        ret[0] = str.s;
        return ret;        
    }
    else {
        if (is_circle) {
            kstring_t str = {0,0,0};
            kputs(s+s2, &str);
            kputsn(s, s1, &str);
            kputs("", &str);
            *n_reads = 1;
            char **ret = malloc(sizeof(void*)*(*n_reads));
            ret[0] = str.s;
            return ret;        
        }
        else {
            char **ret = NULL;          
            if (s2 < l) {
                kstring_t str = {0,0,0};
                kputs(s+s2, &str);
                *n_reads = 1;
                ret = malloc(sizeof(void*)*(*n_reads));
                ret[0] = str.s;                
            }

            if (s1 > 0) {
                kstring_t str = {0,0,0};
                kputsn(s, s1, &str);
                kputs("", &str);
                ret = realloc(ret,sizeof(void*)*(1+*n_reads));
                ret[*n_reads] = str.s;
                *n_reads += 1;
            }

            return ret;
        }
        
    }

    return NULL;
}
static struct hits *check_kmers(struct ref *r, char *s)
{
    khint_t k;
    k = kh_get(hit, r->map, s);
    if (k == kh_end(r->map)) return NULL;
    return kh_val(r->map, k);
}
/*
static char *find_segment_pe(struct ref *ref, struct bseq *seq)
{
    kstring_t str = {0,0,0};
    
}
*/
static int parse_args(int argc, char **argv)
{

    if (argc == 1) return 1;

    const char *thread = NULL;
    const char *seed = NULL;
    const char *tags = NULL;

    int i;   
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-sl") == 0) var = &seed;
        else if (strcmp(a, "-config") == 0) var = &args.config_fname;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tags") == 0) var = &tags;
        else if (strcmp(a, "-pb") == 0) var = &args.phase_tag;
        else if (strcmp(a, "-f") == 0) {
            args.filter_reads = 1;
            continue;
        }
        else if (strcmp(a, "-k") == 0) {
            args.keep_all = 1;
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
        error("Unknown parameter %s", a);
    }

    if (args.config_fname == NULL) error ("Configure file must be set.");
    if (args.input_fname == NULL) error("No input fastq.");

    if (seed) {
        args.seed_length = str2int((char*)seed);
        if (args.seed_length < MIN_SEED_LENGTH || args.seed_length > MAX_SEED_LENGTH)
            error ("Seed length exceed the limit [%d-%d].", MIN_SEED_LENGTH, MAX_SEED_LENGTH);
    }

    if (thread) args.n_thread = str2int((char*)thread);
        
    args.r = config_init(args.config_fname);

    if (args.output_fname)
        args.out = fopen(args.output_fname, "w");
    else
        args.out = stdout;

    if (tags) 
        args.tag_dict = str2tag(tags);
    
    return 0;
}
static void memory_release()
{
    ref_destroy(args.r);    
    fclose(args.out);
}

static struct tag_val *tag_val_create(struct ref *ref)
{
    struct ref_pat *r = ref->ref;
    struct tag_val *v = malloc(sizeof(*v));
    v->tag = dict_init();
    int i;
    for (i = 0; i < r->n; ++i) {
        struct segment *seg = &r->segs[i];
        if (seg->n > 1) {
            int k;
            for (k = 0; k < seg->n; ++k) dict_push(v->tag,seg->s[k].tag);
        }
        else {
            dict_push(v->tag, seg->tag);
        }
    }
    v->n = r->n_tag;
    v->dict = malloc(v->n*sizeof(void *));
    memset(v->dict, 0, sizeof(void*)*v->n);
    
    return v;
}
static void tag_val_destroy(struct tag_val *v)
{
    int i;
    for (i = 0; i < v->n; ++i) {
        if (v->dict[i] != NULL)
            dict_destroy(v->dict[i]);        
    }
    free(v->dict);
    dict_destroy(v->tag);
    free(v);
}
struct ret_dat {
    int m,n;
    struct rd {
        char *name;
        char *seq;
    } *rd;
};
static char *find_segment(struct ref *ref, struct read *r, int n)
{
    struct tag_val *v = tag_val_create(ref);
    struct ret_dat *rd = malloc(sizeof(*rd));
    memset(rd, 0, sizeof(struct ret_dat));
              
    int ir;
    for (ir = 0; ir < n; ++ir) {
        struct read *seq = &r[ir];
        char *s = check_circle(seq->s0);
        kstring_t str = {0,0,0};
        int is_circle = 0;
        if (s) is_circle = 1;        
        else s = strdup(seq->s0);
        int l = strlen(s);
        int length = is_circle ? l + args.seed_length : l - args.seed_length;
        int i;
        for (i = 0; i < length; ++i) {
            str.l = 0;
            int j;
            for (j = 0; j < args.seed_length; ++j) {
                int k = i + j;
                if (k >= l) k -= l;
                kputc(s[k], &str);
            }
            kputs("", &str);
            struct hits *hh = check_kmers(ref, str.s);
            if (hh == NULL) continue;

            for (j = 0; j < hh->n; ++j) {
                struct hit *h = &hh->hit[j];
                int n;
                char **rs =check_pattern(s, i, ref, h, is_circle, v, &n);
                if (n == 0) continue;
                int k;
                for (k = 0; k <n; ++k) {
                    if (rd->n == rd->m) {
                        rd->m = rd->m == 0 ? 2 : rd->m<<1;
                        rd->rd = realloc(rd->rd,rd->m*sizeof(struct rd));
                    }
                    rd->rd[rd->n].name = strdup(seq->name);
                    rd->rd[rd->n].seq = rs[k];
                    rd->n++;
                }
                goto next_read;
            }
        }

        if (rd->n == rd->m) {
            rd->m = rd->m == 0 ? 2 : rd->m<<1;
            rd->rd = realloc(rd->rd,rd->m*sizeof(struct rd));
        }
        rd->rd[rd->n].name = strdup(seq->name);
        rd->rd[rd->n].seq = strdup(seq->s0);
        rd->n++;
        
      next_read:
        free(str.s);
    }

    char *new_name = tagval2name(ref, v);
    
    if (args.keep_all == 0 && new_name == NULL) {
        tag_val_destroy(v);
        return NULL;
    }
    kstring_t out = {0,0,0};
    for (ir = 0; ir < rd->n; ++ir) {
        kputc('>', &out);
        kputs(rd->rd[ir].name, &out);
        if (new_name) kputs(new_name, &out);
        kputc('\n', &out);
        kputs(rd->rd[ir].seq, &out);
        kputc('\n', &out);
        free(rd->rd[ir].name);
        free(rd->rd[ir].seq);
    }
    free(rd->rd);
    free(rd);
    tag_val_destroy(v);
    if (out.l == 0) return NULL;
    return out.s;
}
struct read_block *build_phase_block(struct read_block *r, const char *tag, int *n_block)
{
    if (r->n == 1) {
        *n_block = 1;
        struct read_block *b = read_block_copy(r);
        return b;
    }
    int n = 0, m = 0;
    struct read_block *rb = NULL;
    struct dict *phase_dict = dict_init();
    int i;
    for (i = 0; i < r->n; ++i) {
        struct read *b = &r->b[i];
        char *pv = read_name_pick_tag(b->name, tag);
        if (pv == NULL) { // empty record, treat as
            if (n == m) {
                int old_m = m;
                m = m == 0 ? 2 : m <<1;
                rb = realloc(rb,m*sizeof(struct read_block));
                for (; old_m <m;++old_m)
                    memset(&rb[old_m], 0, sizeof(struct read_block));
            }
            read_block_push(&rb[n], b);
            n++;
        }
        else {
            int id = dict_push(phase_dict, pv);
            if (id >= m) {
                int old_m = m;
                m = id+1;
                rb=realloc(rb, m*sizeof(struct read_block));
                for (; old_m <m;++old_m)
                    memset(&rb[old_m], 0, sizeof(struct read_block));
            }
            read_block_push(&rb[id], b);
            if (n < id) n = id;
        }
    }
    *n_block = n;
    dict_destroy(phase_dict);
    return rb;
}
static void *run_it(void *_p)
{
    struct thread_dat *d = (struct thread_dat*)_p;
    
    kstring_t str = {0,0,0};
    int i; 
    for (i = 0; i < d->n; ++i) {
        struct read_block *rb = &d->rb[i];
        if (args.phase_tag) {
            int n;
            struct read_block *phase_block = build_phase_block(rb, args.phase_tag, &n);
            int j;
            for (j = 0; j < n; ++j) {
                char *s = find_segment(args.r, phase_block[j].b, phase_block[j].n);
                if (s) {
                    kputs(s, &str);
                    free(s);
                }
                int k;
                for (k = 0; k > phase_block[j].n; ++k) {
                    struct read *b = &phase_block[j].b[k];
                    free(b->name);
                    free(b->s0);
                    if (b->q0) free(b->q0);
                    if (b->s1) free(b->s1);
                    if (b->q1) free(b->q1);
                }
            }

        }
        else {
            int j;
            for (j = 0; j < rb->n; ++j) {
                char *s = find_segment(args.r, &rb->b[j], 1);
                if (s) {
                    kputs(s, &str);
                    free(s);
                }
            }
        }
    }

    thread_dat_destroy(d);
    if (str.l == 0) return NULL;
    return str.s;
}
static void write_out(void *_d)
{
    char *s = (char*)_d;
    if (s != NULL && strlen(s) > 0) {
        fprintf(args.out, "%s", s);
        free(s);
    }
}

int fastq_segment(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv) == 1) return usage();

    FILE *fp_in = fopen(args.input_fname, "r");
    if (fp_in == NULL) error("%s : %s.", args.input_fname, strerror(errno));
   
    FILE *fp_out = stdout;
    if (args.output_fname) {
        fp_out = fopen(args.output_fname, "w");
        if (fp_out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    }
    
    LOG_print("Build reference finished.");

    if (args.n_thread == 1) {
        for (;;) {
            struct thread_dat *dat = read_thread_dat(fp_in, args.tag_dict);
            if (dat == NULL) break;
            char *ret = (char*)run_it(dat);
            if (ret == NULL) continue;
            fputs(ret, fp_out);
            free(ret);
        }
    }
    else {
        hts_tpool *p = hts_tpool_init(args.n_thread);
        hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
        hts_tpool_result *r;

        for (;;) {
            struct thread_dat *dat = read_thread_dat(fp_in, args.tag_dict);
            if (dat == NULL) break;
            
            int block;
            do {
                block = hts_tpool_dispatch2(p, q, run_it, dat, 1);
                if ((r = hts_tpool_next_result(q))) {
                    char *d = (char*) hts_tpool_result_data(r);
                    if (d) {
                        fputs(d, fp_out);
                        free(d);
                    }
                    hts_tpool_delete_result(r, 0);
                }
            } while(block == -1);
        }
        hts_tpool_process_flush(q);

        while ((r = hts_tpool_next_result(q))) {
            char *d = (char*) hts_tpool_result_data(r);
            if (d) {
                fputs(d, fp_out);
                free(d);
            }
            hts_tpool_delete_result(r, 0);
        }
        hts_tpool_process_destroy(q);
        hts_tpool_destroy(p);
    }

    fclose(fp_in);
    fclose(fp_out);
    memory_release();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
/*
int main(int argc, char **argv)
{
    return check_segment2(argc, argv);
}
*/
