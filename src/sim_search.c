#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "sim_search.h"

#define BASE_TERM 0x0
#define BASE_A  0x1
#define BASE_C  0x2
#define BASE_G  0x3
#define BASE_T  0x4
#define BASE_N  0x5

static int kmer_min = 5;
static int kmer_max = 21;

static int kmer_size = 5;

typedef struct ss_idx {
    int *idx;
    int n,m;
} sidx_t;

KHASH_MAP_INIT_INT64(ss64, int)
KHASH_MAP_INIT_INT(ss32, sidx_t)

typedef kh_ss64_t hash64_t;
typedef kh_ss32_t hash32_t;

struct similarity_search_aux {
    hash64_t *d0;
    hash32_t *d1; 
    uint64_t *cs; // compact sequence    
    int n, m;
};

uint8_t encode_base(char c)
{
    switch (c) {
        case 'A':
        case 'a':
            return BASE_A;
        case 'C':
        case 'c':
            return BASE_C;
        case 'G':
        case 'g':
            return BASE_G;
        case 'T':
        case 't':
            return BASE_T;
        case 'N':
        case 'n':
            return BASE_N;
        case '\0':
            return 0x0;
        default:
            error("Try to encode a non DNA sequence ? %c", c);
    }
}
// not safe for Ns
uint32_t enc32(char *s, int l)
{
    int len = strlen(s);
    if (len < l) error("Try to encode a truncated sequence ?");
    if (len > 16) error("Only support to encode sequence not longer than 16nt.");
    uint64_t q = 0;
    int i;
    for (i = 0; i < l; ++i)
        q = q<<3 | (encode_base(s[i]) & 0x7);
    return q;

}
// not safe for Ns
uint64_t enc64(char *s)
{
    int l = strlen(s);
    if (l > 32) error("Only support to encode sequence not longer than 32nt.");
    uint32_t q = 0;
    int i;
    for (i = 0; i < l; ++i)
        q = q<<3 | (encode_base(s[i]) & 0x7);
    return q;
}

char *decode64(uint64_t q)
{
    char c[23];
    memset(c, 0, 23);
    int i = 21;
    for (;;) {
        uint8_t x = q & 0x7;
        if (x == 0x1) c[i] = 'A'; // kputc('A', &str);
        else if (x == 0x2) c[i] = 'C'; // kputc('C', &str);
        else if (x == 0x3) c[i] = 'G'; // kputc('G', &str);
        else if (x == 0x4) c[i] = 'T'; // kputc('T', &str);
        else if (x == 0x5) c[i] = 'N';
        else break;
        i--;
        q = q>>3;
    }
    kstring_t str = {0,0,0};
    kputs(c+i+1, &str);
    return str.s;
}
char *decode32(uint32_t q)
{
    char c[11];
    memset(c, 0, 11);
    int i = 9;
    for (;;) {
        uint8_t x = q & 0x7;
        if (x == 0x1) c[i] = 'A'; // kputc('A', &str);
        else if (x == 0x2) c[i] = 'C'; //kputc('C', &str);
        else if (x == 0x3) c[i] = 'G'; //kputc('G', &str);
        else if (x == 0x4) c[i] = 'T'; //kputc('T', &str);
        else if (x == 0x5) c[i] = 'N';
        else break;
        i--;
        q = q>>3;
    }
    kstring_t str = {0,0,0};
    kputs(c+i+1, &str);
    return str.s;
}
static int check_Ns(char *s, int l)
{
    if (l == 0) l = strlen(s);
    int i;
    for (i = 0; i < l; ++i)
        if (s[i] != 'A' && s[i] != 'a' && s[i] != 'C' && s[i] != 'c'
            && s[i] != 'G' && s[i] != 'g' && s[i] != 'T' && s[i] != 't') return 1;

    return 0;
}
ss_t *ss_init()
{
    ss_t *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(ss_t));
    s->d0 = kh_init(ss64);
    s->d1 = kh_init(ss32);
    return s;
}

void ss_destroy(ss_t *S)
{
    kh_destroy(ss64, S->d0);
    khint_t k;
    for (k = kh_begin(S->d1); k != kh_end(S->d1); ++k) {
        if (kh_exist(S->d1, k)) {
            sidx_t *idx = &kh_val(S->d1, k);
            if (idx &&idx->idx) free(idx->idx);
        }
    }
    kh_destroy(ss32, S->d1);
    free(S->cs);
    free(S);
}
static char *print_bits64(uint64_t x)
{
    kstring_t str={0,0,0};
    int i;
    for(i = 0;i<64;++i){
        kputc(x&((uint64_t)1<<63) ? '1': '0', &str);
        x = x<<1;
    }
    return str.s;
}
static char *print_bits32(uint32_t x)
{
    kstring_t str={0,0,0};
    int i;
    for(i = 0;i<32;++i){
        kputc(x&((uint32_t)1<<31) ? '1': '0', &str);
        x = x<<1;
    }
    return str.s;
}
static void build_kmers(ss_t *S, uint64_t q, int idx)
{
    int offset = 3 *kmer_size;
    uint32_t mask = ~(0x1<<offset);
    mask = mask<<(32-offset)>>(32-offset);
    for (;;) {
        if (q>>(offset-3) == 0) break;
        uint32_t x = q & mask;
        q=q>>3;
        khint_t k = kh_get(ss32, S->d1, x);
        if (k != kh_end(S->d1)) {
            struct ss_idx *si = &kh_val(S->d1, k);
            if (si->m == si->n) {
                si->m = si->m<<1;
                si->idx = realloc(si->idx, si->m*sizeof(int));
            }
            si->idx[si->n++] = idx;
        }
        else {
            int ret;
            k = kh_put(ss32, S->d1, x, &ret);
            struct ss_idx *si = &kh_val(S->d1, k);
            memset(si, 0, sizeof(struct ss_idx));
            si->m = 2;
            si->idx = realloc(si->idx, sizeof(int)*si->m);
            si->idx[si->n++] = idx;
        }
 
    }
}
int ss_push(ss_t *S, char *seq)
{
    int N = check_Ns(seq, 0);
    if (N) error("Try to push sequence %s contain Ns.", seq);
    uint64_t q = enc64(seq);
    khint_t k = kh_get(ss64, S->d0, q);
    if (k != kh_end(S->d0)) return 1;
    if (S->n == S->m) {
        S->m = S->m == 0 ? 1024 : S->m<<1;
        S->cs = realloc(S->cs, S->m*sizeof(uint64_t));
    }
    S->cs[S->n] = q;
    int ret;
    k = kh_put(ss64, S->d0, S->cs[S->n], &ret);
    kh_val(S->d0, k) = S->n;
    build_kmers(S, q, S->n);
    S->n++;
    return 0;
}

struct element {
    int ele;
    int cnt;
};
typedef struct set {
    struct element *ele;
    int n, m;
} set_t;

set_t *set_init()
{
    set_t *set = malloc(sizeof(*set));
    memset(set, 0, sizeof(*set));
    return set;
}
void set_destory(set_t *set)
{
    if (set->m) free(set->ele);
    free(set);
}
static void set_push_core(int ele, set_t *set)
{
    int i;
    for (i = 0; i < set->n; ++i)
        if (set->ele[i].ele == ele) {
            set->ele[i].cnt++;
            return;
        }

    if (set->n == set->m) {
        set->m = set->m == 0 ? 2 : set->m<<1;
        set->ele = realloc(set->ele, sizeof(struct element)*set->m);
    }
    
    set->ele[set->n].ele = ele;
    set->ele[set->n].cnt = 1;
    set->n++;
}
void set_push(int *ele, int n, set_t *set)
{
    assert(n > 0);

    int i;
    for (i = 0; i < n; ++i) set_push_core(ele[i], set);    
}
int cmpfunc (const void * a, const void * b)
{
    return (*(struct element*)b).cnt - (*(struct element*)a).cnt;
}

int set_top_2(set_t *set)
{
    if (set->n <= 2) return set->n;

    qsort(set->ele, set->n,  sizeof(struct element),  cmpfunc);

    int max = set->ele[1].cnt;
    int i;
    for (i = 2; max <= set->ele[i].cnt && i < set->n; ++i) {}
    return i;
}

extern int levenshtein(char *a, char *b, int l);
int levnshn_dist_calc(uint64_t a, uint64_t b)
{
    char *s1 = decode64(a);
    char *s2 = decode64(b);
    int l = strlen(s1);
    int dist = levenshtein(s1, s2, l);
    free(s1);
    free(s2);
    return dist;
}

char *ss_query(ss_t *S, char *seq, int e, int *exact)
{
    *exact = 1; // exactly match
    int l = strlen(seq);
    if (l > kmer_max) error("Sequence is too long. %s", seq);
    uint64_t q = enc64(seq);
    khint_t k = kh_get(ss64, S->d0, q);
    
    if (k != kh_end(S->d0)) return decode64(q);

    *exact = 0;
    int i;
    set_t *set = set_init();
    
    for (i = 0; i < l - kmer_size; ++i) {
        int j = check_Ns(seq+i, kmer_size);
        if (j) {
            i+=j;
            continue;
        }
        uint32_t q0 = enc32(seq+i, kmer_size);
        // char *s0 = decode32(q0);

        k = kh_get(ss32, S->d1, q0);
        if (k == kh_end(S->d1)) continue;
        
        struct ss_idx *idx = &kh_val(S->d1, k);
        set_push(idx->idx, idx->n, set);
    }

    int n = set_top_2(set);
    int hit = -1;
    for (i = 0; i < n; ++i) {
        int dist = levnshn_dist_calc(S->cs[set->ele[i].ele], q);
        if (dist <= e) {
            if (hit != -1) goto multi_hits;
            hit = set->ele[i].ele;
        }
    }
    
    set_destory(set);

    if (hit == -1) return NULL;
    return decode64(S->cs[hit]);

  multi_hits:
    set_destory(set);
    return NULL;
}
/*
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

*/
#ifdef SS_MAIN
int main()
{
    ss_t *S = ss_init();
    ss_push(S,"TAACACGCAA");
    ss_push(S,"TAACAGCCAA");
    char *s1 = ss_query(S, "TAACAGCCAA", 1);
    char *s2 = ss_query(S, "TAACAGCCNA", 2);
    fprintf(stderr, "%s\t%s\n", s1, s2);
    return 0;
}

#endif
