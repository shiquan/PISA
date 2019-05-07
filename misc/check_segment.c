#include "utils.h"
#include "json_config.h"
#include "fastq.h"
#include "kson.h"
#include "number.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/thread_pool.h"

#define MAX_PATTERN_LENGTH 100
#define MIN_SEED_LENGTH    6
#define MAX_SEED_LENGTH    20
#define MAX_MISMATCH       0

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

struct segment {
    int l;
    int dist;
    int n;
    char **wl; // white list
    kh_str_t *hash; // white list hash;    
};

struct ref_pat {
    int idx[MAX_PATTERN_LENGTH]; // fast access the index of segment
    int pos[MAX_PATTERN_LENGTH]; // the position of each segment, start from 0
    char *seq; // ACGTN, N for tag
    int l_seq;

    int n;
    struct segment *segs;
};
struct ref {
    struct ref_pat *ref;
    struct ref_pat *rev;
    int n, m;
    char **kmers;
    kh_hit_t *map;
};
struct ref_pat *ref_pat_alloc()
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
        if (s->n > 0) {
            int j;
            for (j = 0; j < s->n; ++j) free(s->wl[j]);
            free(s->wl);
            kh_destroy(str,s->hash);
        }
        // if (s->str) free(s->str);
    }
    free(r->seq);
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
    free(r);
}
static char *rev_seq(char *s, int l)
{
    char *r = malloc(l*sizeof(char));
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
            default:
                error("Try to reverse a non DNA sequence? \"%s\"", s);
        }
    }
    return r;
}
static struct ref_pat *build_rev_pat(struct ref_pat *r)
{
    struct ref_pat *v = ref_pat_alloc();
    int i;
    v->n = r->n;
    v->l_seq = r->l_seq;
    v->seq = malloc(v->l_seq*sizeof(char));
    
    int p = 0;
    
    for (i = 0; i < v->l_seq; ++i) {
        switch(r->seq[v->l_seq-i-1]) {
            case 'A':
                v->seq[i] = 'T';
                break;
            case 'C':
                v->seq[i] = 'G';
                break;
            case 'G':
                v->seq[i] = 'C';
                break;
            case 'T':
                v->seq[i] = 'A';
                break;
            default:
                v->seq[i] = 'N';
                break;
        }

        v->idx[i] = r->seq[v->l_seq-i-1];

        if (v->idx[i] == 0) p = 0;
        else v->pos[i] = p++;
    }

    v->segs = malloc(v->n*sizeof(struct segment));
    for (i = 0; i < v->n; ++i) {
        struct segment *ori = &r->segs[i];
        struct segment *copy = &v->segs[i];
        copy->l = ori->l;
        copy->dist = ori->dist;
        copy->n = ori->n;
        if (copy->n) {
            copy->hash = kh_init(str);
            copy->wl = malloc(sizeof(char*)*copy->n);
            int j, k;
            khint_t ik;
            int ret;
            for (j = 0; j < copy->n; ++j) {
                copy->wl[j] = malloc(sizeof(char)*copy->l);
                char *c = ori->wl[j];
                for (k = 0; k < copy->l; ++k) {
                    switch(c[copy->l -k -1]) {
                        case 'A':
                            copy->wl[j][k] = 'T';
                            break;
                        case 'C':
                            copy->wl[j][k] = 'G';
                            break;
                        case 'G':
                            copy->wl[j][k] = 'C';
                            break;
                        case 'T':
                            copy->wl[j][k] = 'A';
                            break;
                        default:
                            error("Unknown character at white list, %c", c[copy->l-k-1]);
                    }
                }
                ik = kh_put(str, copy->hash, copy->wl[j], &ret);
                kh_val(copy->hash, ik) = j;
            }
        }
    }
    return v;
}
static struct args {
    const char *r1_fname;
    const char *r2_fname;
    const char *bam_fname;
    const char *barcode_tag;
    const char *config_fname;
    const char *output_fname;

    int seed_length;

    int smart_pairing;
    struct ref *r;

    struct fastq_handler *fastq;
    
    FILE *report_fp;

    uint64_t found;
    uint64_t partly_found;
    
    int n_thread;
} args = {
    .r1_fname = NULL,
    .r2_fname = NULL,
    .bam_fname = NULL,
    .barcode_tag = NULL,
    .config_fname = NULL,
    .output_fname = NULL,

    .seed_length = 6,

    .smart_pairing = 0,
    .r = NULL,

    .fastq = NULL,
    
    .report_fp = NULL,

    .found = 0,
    .partly_found = 0,

    .n_thread = 1,
};


static struct ref_pat *config_init(const char *fn)
{
    char *config_str = json_config_open(fn);
    if (config_str == NULL) error("Empty configure file.");
    kson_t *json = kson_parse(config_str);
    free(config_str);

    char *pat = NULL;

    struct ref_pat *ref = ref_pat_alloc();

    kh_str_t *hash = kh_init(str);
    char **tags = NULL;
    
    do {
        const kson_node_t *root = json->root;
        if (root == NULL) error("Format error. Root node is emtpy.");
        int i;
        for (i = 0; i < root->n; ++i) {
            const kson_node_t *node = kson_by_index(root, i);
            if (node == NULL) continue;
            if (node->key == NULL) error("Format error. Node key is empty.");
            if (strcmp(node->key, "pattern") == 0) {
                pat = strdup(node->v.str);
            }
            else if (strcmp(node->key, "segments") == 0) {

                if (node->type != KSON_TYPE_BRACKET) error("Format error.\"segments\":[{},{}]");

                ref->n = node->n;
                if (ref->n == 0) error("No segment records in the config file.");
                ref->segs = malloc(ref->n*sizeof(struct segment));
                tags = malloc(ref->n*sizeof(char*));
                
                int j;
                for (j =0; j <node->n; ++j) {
                    const kson_node_t *n = kson_by_index(node, j);
                    if (n == NULL) error("Segment record is empty.");
                    if (n->type!= KSON_TYPE_BRACE) error("Format error.\"segments\":[{},{}]");
                    struct segment *seg = &ref->segs[j];
                    memset(seg, 0, sizeof(struct segment));

                    int k;
                    for (k = 0; k < n->n; ++k) {
                        const kson_node_t *n1 = kson_by_index(n, k);
                        if (n1 == NULL) error("Segment record is empty.");
                        
                        if (strcmp(n1->key,"tag") == 0) {
                            tags[j] = strdup(n1->v.str);
                        }
                        else if (strcmp(n1->key, "length") == 0) {
                            seg->l = str2int(n1->v.str);
                        }
                        else if (strcmp(n1->key, "distance") == 0) {
                            seg->dist = str2int(n1->v.str);
                        }
                        else if (strcmp(n1->key, "white list") == 0) {
                            if (node->type != KSON_TYPE_BRACKET) error("Format error.\"white list\":[]");
                            seg->n = n1->n;
                            seg->wl = malloc(seg->n*sizeof(char*));
                            seg->hash = kh_init(str);
                            int m;
                            for (m = 0; m <n1->n; ++m) {
                                const kson_node_t *n2 = kson_by_index(n1,m);
                                if (n2->v.str) error ("Empty record at white list.");
                                seg->wl[m] = strdup(n2->v.str);
                                khint_t k;
                                int ret;
                                k = kh_put(str, seg->hash, seg->wl[m], &ret);
                                if (!ret) error ("Duplicate record at white list. %s", seg->wl[m]);
                                kh_val(seg->hash, k) = m;                                
                            }
                        }
                    }
                    // check the value
                    if (seg->l == 0) error("\"length\" must be set for each segment.");
                    if (tags[j] == NULL) error("\"tag\" must be set for each segment.");
                    khint_t ik;
                    int ret;
                    ik = kh_put(str, hash, tags[j], &ret);
                    if (!ret) error("Tag %s is duplicated at configure file.", tags[j]);
                    kh_val(hash, ik) =j;
                }
            }
            else error("Unknown key, %s", node->key);
        }
    } while(0);

    kson_destroy(json);
    
    if (pat == NULL) error("\"pattern\" must be set in the configure file.");
    
    int l = strlen(pat);
    if (l > MAX_PATTERN_LENGTH) error("pattern sequence is too long. Max length is %d.", MAX_PATTERN_LENGTH);

#define _fill_id(start, length, fill) do {\
        int _i, _j;\
        for (_i = start, _j = 0; _j < length; ++_i,++_j) ref->idx[_i] = fill; \
    } while(0)
    
    do {

        kstring_t str = {0,0,0};
        kstring_t seq = {0,0,0};
        int is_tag = 0;
        int i;
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
                        is_tag = 1;
                    } while(0);
                    break;
                    
                case ')':
                    do {
                        if (is_tag == 0) error("Unknown format.");
                        if (str.l != 2) error("Tag name is too long. Only support 2 character, such 'AS','CB' etc.");
                        
                        khint_t k;
                        k = kh_get(str, hash, str.s);
                        if (k == kh_end(hash)) error("No tag %s found at configure.", str.s);
                        int id = kh_val(hash, k);
                        struct segment *s = &ref->segs[id];
                        if (seq.l + s->l > MAX_PATTERN_LENGTH) error("Pattern sequence is too long. Max is %d.", MAX_PATTERN_LENGTH);
                        
                        _fill_id(seq.l, s->l, id);

                        int _i;
                        for (_i = 0; _i < s->l; ++_i) kputc('N', &seq);
                        
                        str.l = 0;
                        is_tag = 0;
                    } while(0);
                    break;
                    
                case 'a': case 'A':
                case 'c': case 'C':                
                case 'g': case 'G':
                case 't': case 'T':
                    kputc(pat[i], &str);
                    break;
                        
                case 'b': case 'B':                
                case 'd': case 'D':
                case 'e': case 'E':
                case 'f': case 'F':
                case 'h': case 'H':
                case 'i': case 'I':
                case 'j': case 'J':
                case 'k': case 'K':
                case 'l': case 'L':
                case 'm': case 'M':
                case 'n': case 'N':
                case 'o': case 'O':
                case 'p': case 'P':
                case 'q': case 'Q':
                case 'r': case 'R':
                case 's': case 'S':
                case 'u': case 'U':
                case 'v': case 'V':
                case 'w': case 'W':
                case 'x': case 'X':
                case 'y': case 'Y':
                case 'z': case 'Z':
                    do {
                        if (is_tag == 0) error("Unsupport DNA sequence in the pattern string. '%c'", pat[i]);
                        kputc(pat[i], &str);
                        if (str.l > 2) error("Tag name is too long. Only support 2 character, such 'AS','CB' etc.");
                    } while(0);
                    break;
                    
                default:
                    error("Unsupport character '%c' in the pattern string.", pat[i]);
                    break;
            }
        }

#undef _fill_id

        if (is_tag == 1) error("Unknown format.");
        if (str.l) kputs(str.s, &seq);
        free(str.s);

        ref->seq = seq.s;
        ref->l_seq = seq.l;

        for (i = 0; i < ref->l_seq; ++i) {
            switch(ref->seq[i]) {
                case 'a':
                    ref->seq[i] = 'A';
                    break;

                case 'c':
                    ref->seq[i] = 'C';
                    break;

                case 'g':
                    ref->seq[i] = 'G';
                    break;

                case 't':
                    ref->seq[i] = 'T';
                    break;

                default:
                    break;
            }
        }
                        
    } while (0);

    int i;
    fputs("#RNAME",args.report_fp);
    for (i = 0; i < ref->n; ++i) {
        fprintf(args.report_fp, "\t%s", tags[i]);
        free(tags[i]);
    }
    fputc('\n', args.report_fp);
    free(tags);
    kh_destroy(str, hash);
        
    free(pat);
    
    return ref;    
}
static void build_kmers_ref(struct ref *r)
{
    char *ref[2];
    ref[0] = r->ref->seq;
    ref[1] = r->rev->seq;

    r->map = kh_init(hit);
    
    int strand = 0;
    for ( ;strand < 2; ++strand) {
        int i, j;
        khint_t k;
        int ret;

        for (i = 0; i < r->ref->l_seq-args.seed_length+1;) {
            char *s = ref[strand];
            if (s[i] == 'N') { i++; continue; }
            for (j = 0; j < args.seed_length && i+j < r->ref->l_seq; j++)
                if (s[i+j] == 'N') break;
            if (j == args.seed_length) {
                kstring_t str = {0,0,0};
                kputsn(s+i, args.seed_length, &str);
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
                }
                else {
                    struct hits *h = kh_val(r->map, k);
                    h->hit = realloc(h->hit, (h->n+1)*sizeof(struct hit));
                    struct hit *hh = &h->hit[h->n++];
                    hh->strand = strand;
                    hh->loc = i;
                    free(str.s); // free key
                }
            }
            i++;
        }
    }
}
static int usage()
{
    fprintf(stderr, "* Pick pre-designed segments from reads, export CB code and mark code to a table.\n");
    fprintf(stderr, "PickSegFromReads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "-1      [fastq]   Read 1 or smart pairing reads.\n");
    fprintf(stderr, "-2      [fastq]   Read 2.\n");
    // fprintf(stderr, "-bam    [bam]     Input bam file, sorted by read name. Usually be unmapped reads. Conflict with -1 and -2.\n");
    fprintf(stderr, "-p                Input is smart pairing. Conflict with -2.\n");
    fprintf(stderr, "-cb     [CB]      Cell barcode tag in the read name, processed by SingleCellTools parse.\n");
    fprintf(stderr, "-config [json]    Configure file.\n");    
    fprintf(stderr, "-o      [tsv]     Output table.\n");
    fprintf(stderr, "-sl     [INT]     Seed length for mapping consensus sequence.\n");
    return 1;
}

/*
static kstring_t *build_fragment_from_PE()
{
}
*/

#define FQ_FLG_UNK 0
#define FQ_FLG_FUD 1
#define FQ_FLG_PAR 2

struct pos_p {
    int id;
    int loc;
};
// n is the count of hits, 0 for none, -1 for bad matched (mismatches exceed edit distance)
static struct pos_p *check_pattern_forward(char *s, const int start, const struct ref_pat *r, const struct hit *h, int *n)
{
    int i = h->loc+args.seed_length; // for ref
    int c = 0;
    struct pos_p *p = malloc(sizeof(struct pos_p));
    int loc = start + args.seed_length; // for sequence
    int edit_dist = 0;
    int len = strlen(s);
    
    *n = 0;
    for (; i < r->l_seq && loc < len;) {
        if (r->seq[i] == 'N') {
            p[c].loc = loc;
            p[c].id = r->idx[i];
            c++;
            p = realloc(p, sizeof(struct pos_p)*c);
            for(;r->seq[i]=='N'; ++i, ++loc);
            continue;
        }
        
        if (s[loc] != r->seq[i]) edit_dist++;
        if (edit_dist > MAX_MISMATCH) {
            *n = -1;
            free(p);
            return NULL; 
        }        
        loc++;
        i++;
    }
    if (c == 0) {
        free(p);
        return NULL;
    }
    *n = c;    
    return p;
}
static struct pos_p *check_pattern_backward(char *s, const int start, const struct ref_pat *r, const struct hit *h, int *n)
{
    int i = h->loc;
    int c = 0;
    struct pos_p *p = malloc(sizeof(struct pos_p));
    int loc = start;
    int edit_dist = 0;
    
    *n = 0;
    for (; i >= 0 && loc >= 0;) {
        if (r->seq[i] == 'N') {
            for(; i>=0 && r->seq[i]=='N'; --i, --loc);
            if (loc >= 0) {
                p[c].loc = loc+1;
                p[c].id = r->idx[i+1];            
                c++;
                p = realloc(p, sizeof(struct pos_p)*c);
            }
            continue;
        }
        
        if (s[loc] != r->seq[i]) edit_dist++;
        if (edit_dist > MAX_MISMATCH) {
            *n = -1;
            free(p);
            return NULL; 
        }        
        loc--;
        i--;
    }
    *n = c;
    return p;

}
static int check_segment(struct segment *G, char *s)
{
    khint_t k;
    k = kh_get(str, G->hash,s);
    if (k == kh_end(G->hash)) return -1;
    return kh_val(G->hash,k);
}

static char **check_pattern(char *s, int start, struct ref_pat *r, struct hit *h, int *partly_found)
{
    int l1, l2;
    struct pos_p *p1 = check_pattern_forward(s, start, r, h, &l1);
    if (l1 == -1) return NULL;
    
    struct pos_p *p2 = check_pattern_backward(s, start, r, h, &l2);
    if (l2 == -1) return NULL;
    
    if (l1 + l2 < r->n) {
        if (p1) free(p1);
        if (p2) free(p2);
        *partly_found = 1;
        return NULL;
    }

    int len = strlen(s);
    
    char **fetch = malloc(sizeof(char*)*r->n);
    memset(fetch, 0, sizeof(char*)*r->n);
    
    int i;
    for (i = 0; i < r->n; ++i) {
        kstring_t str = {0,0,0};
        struct pos_p *p = NULL;
        if (i >= l1) p = &p2[i-l1];
        else p = &p1[i];            
        struct segment *seg = &r->segs[p->id];
        if (p->loc + seg->l > len) { //partly found
            int j;
            for (j = 0; j < r->n; ++j)
                if (fetch[j]) free(fetch);
            free(fetch);
            return NULL;
        }
        
        kputsn(s+p->loc, seg->l, &str); kputs("", &str);
        
        if (seg->n > 0) {
            int e = check_segment(seg, str.s);
            if (e < 0) { // unmatched
                int j;
                for (j = 0; j < r->n; ++j)
                    if (fetch[j]) free(fetch);
                free(fetch);
                free(str.s);
                return NULL;
            }
            str.l = 0;
            kputs(seg->wl[e], &str);
        }

        if (h->strand == 1) { // reverse sequence for minus strand
            char *r = rev_seq(str.s, str.l);
            str.l = 0;
            kputs(r, &str);
            free(r);
        }
        
        // Concat segment with same tag name into one
        if (fetch[p->id] != NULL) {
            kputs(fetch[p->id], &str);
            free(fetch[p->id]);
        }
        
        fetch[p->id] = str.s;        
    }
    
    return fetch;
}
static struct hits *check_kmers(struct ref *r, char *s)
{
    khint_t k;
    k = kh_get(hit, r->map, s);
    if (k == kh_end(r->map)) return NULL;
    return kh_val(r->map, k);
}
// @l      Length of string
static char **find_segment_core(struct ref *r, char *s, int l, int *partly_found)
{
    int i;
    kstring_t str ={0,0,0};
    int length = l - args.seed_length;
    
    for (i = 0; i < length; ++i) {
        kputsn(s+i, args.seed_length, &str);
        struct hits *hh = check_kmers(r, str.s);
        if (hh) {
            int j;
            for (j = 0; j < hh->n; ++j) {
                struct hit *h = &hh->hit[j];
                char **fetch = check_pattern(s, i, h->strand == 0 ? r->ref : r->rev, h, partly_found);
                if (fetch) {
                    free(str.s);
                    return fetch;
                }
            }
        }
        str.l = 0;
    }
    free(str.s);
    return NULL;
}
static void find_segment(struct ref *ref, struct bseq *seq)
{
    int partly_found = 0;
    
    char **fetch =  NULL;
    fetch = find_segment_core(ref, seq->s0, seq->l0, &partly_found);
    if (fetch == NULL)
        fetch = find_segment_core(ref, seq->s1, seq->l1, &partly_found);
    if (fetch) {
        kstring_t str = {0,0,0};
        int i;
        for (i = 0; i < ref->ref->n; ++i) {
            if (i) kputc('\t', &str);
            if (fetch[i]) {
                kputs(fetch[i], &str);
                free(fetch[i]);
            }
            else {
                kputc('.', &str);
            }
        }
        free(fetch);
        seq->data = (void *)str.s;
        seq->flag = FQ_FLG_FUD;
    }
    else if (partly_found) seq->flag = FQ_FLG_PAR;
}

static void *run_it(void *_p)
{
    struct bseq_pool *p = (struct bseq_pool*)_p;
    int i;
    for (i = 0; i < p->n; ++i) 
        find_segment(args.r, &p->s[i]);
    
    return p;
}
static void write_out(void *_d)
{
    struct bseq_pool *p = (struct bseq_pool*)_d;
    struct args *opts = (struct args*)p->opts;
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        if (b->flag == FQ_FLG_FUD) {
            fprintf(opts->report_fp, "%s\t%s\n", b->n0, (char*)b->data);
            free(b->data);
            opts->found++;
        }
        else if (b->flag == FQ_FLG_PAR) {
            opts->partly_found++;
        }
    }
    bseq_pool_destroy(p);
}

static int parse_args(int argc, char **argv)
{

    if (argc == 1) return 1;

    int i;
    const char *thread = NULL;
    const char *seed = NULL;

    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return usage();

        if (strcmp(a, "-1") == 0) var = &args.r1_fname;
        else if (strcmp(a, "-2") == 0) var = &args.r2_fname;
        else if (strcmp(a, "-cb") == 0) var = &args.barcode_tag;
        else if (strcmp(a, "-sl") == 0) var = &seed;
        else if (strcmp(a, "-config") == 0) var = &args.config_fname;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-p") == 0) {
            args.smart_pairing = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }                                           
    }

    if (args.config_fname == NULL) error ("Configure file must be set.");
    if (args.output_fname == NULL) error ("Output file must be set.");
    if (args.r1_fname == NULL) error ("Input fastq must be specified with -1.");
    args.fastq = fastq_handler_init(args.r1_fname, args.r2_fname, args.smart_pairing, 100000);
    CHECK_EMPTY (args.fastq, "Failed to init fastq file.");

    args.report_fp = fopen(args.output_fname, "w");
    CHECK_EMPTY(args.report_fp, "%s : %s.", args.output_fname, strerror(errno));

    if (seed) {
        args.seed_length = str2int((char*)seed);
        if (args.seed_length < MIN_SEED_LENGTH || args.seed_length > MAX_SEED_LENGTH)
            error ("Seed length exceed the limit [%d-%d].", MIN_SEED_LENGTH, MAX_SEED_LENGTH);
    }

    if (thread) args.n_thread = str2int((char*)thread);
        
    args.r = malloc(sizeof(struct ref));
    struct ref *r = args.r;
    memset(r, 0, sizeof(*r));
    r->ref = config_init(args.config_fname);
    r->rev = build_rev_pat(r->ref);
    build_kmers_ref(r);
    
    return 0;
}
static void memory_release()
{
    fastq_handler_destory(args.fastq);
    ref_destroy(args.r);    
    fclose(args.report_fp);    
}
int main(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv) == 1) return usage();
    LOG_print("Build reference finished.");
    
    for (;;) {
        struct bseq_pool *p = bseq_read(args.fastq, &args);
        if (p == NULL) break;
        run_it(p);
        write_out(p);
    }

    memory_release();
    
    return 0;
}
