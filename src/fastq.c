#include "fastq.h"
#include "number.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

int check_name(char *s1, char *s2)
{
    int l1 = strlen(s1);
    int l2 = strlen(s2);
    if ( l1 != l2 ) return 1;

    size_t n;
    for(n = 0; n < l1; ++n, ++s1, ++s2) {
        if (*s1 != *s2) return 1;
    }
    return 0;
}

struct bseq_pool *bseq_pool_init()
{
    struct bseq_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    return p;
}
void bseq_pool_clean(struct bseq_pool *p)
{
    int i;
    for ( i = 0; i < p->n; ++i ) {
        struct bseq *b = &p->s[i];
        if (b->l0) {
            free(b->n0);
            free(b->s0);
            if (b->q0) free(b->q0);
        }
        if (b->l1) {
            // free(b->n1);
            free(b->s1);
            if (b->q1) free(b->q1);
        }
    }
    if (p->m > 0) free(p->s);
}
void bseq_pool_destroy(struct bseq_pool *p)
{
    bseq_pool_clean(p);
    free(p);
}

void trim_read_tail(char *s, int l)
{
    if ( l > 2 && s[l-2] == '/' ) s[l-2] = '\0';    
}

static struct bseq_pool *bseq_read_smart(kseq_t *ks, int chunk_size)
{
    struct bseq_pool *p = bseq_pool_init();
    int size = 0;
    do {
        if ( kseq_read(ks) < 0 ) break;
        struct bseq *s;
        if (p->n >= p->m ) {
            p->m = p->m ? p->m*2 : 256;
            p->s = realloc(p->s, p->m*sizeof(struct bseq));
        }
        s = &p->s[p->n];
        memset(s, 0, sizeof(*s));
        trim_read_tail(ks->name.s, ks->name.l);
        s->n0 = strdup(ks->name.s);
        s->s0 = strdup(ks->seq.s);
        s->q0 = ks->qual.l? strdup(ks->qual.s) : 0;
        s->l0 = ks->seq.l;
        size += s->l0;
        if ( kseq_read(ks) < 0 ) error("Truncated input.");

        trim_read_tail(ks->name.s, ks->name.l);
        // s->n1 = strdup(ks->name.s);
        if ( check_name(s->n0, ks->name.s) ) error("Inconsistance paired read names. %s vs %s.", s->n0, ks->name.s);
        s->s1 = strdup(ks->seq.s);
        s->q1 = ks->qual.l? strdup(ks->qual.s) : 0;
        s->l1 = ks->seq.l;
        size += s->l1;
        p->n++;
        
        if ( size >= chunk_size ) break;
    }
    while (1);
    if ( p->n == 0 ) {
        bseq_pool_destroy(p);
        return NULL;
   }
    return p;
}
static struct bseq_pool *bseq_read_core(kseq_t *k1, kseq_t *k2, int chunk_size, int pe)
{
    struct bseq_pool *p = bseq_pool_init();
    if ( pe == 0 ) {
        do {
            if ( kseq_read(k1) < 0 ) break;
            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct bseq));
            }
            struct bseq *s = &p->s[p->n];
            memset(s, 0, sizeof(*s));
            trim_read_tail(k1->name.s, k1->name.l);
            s->n0 = strdup(k1->name.s);
            s->s0 = strdup(k1->seq.s);
            s->q0 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l0 = k1->seq.l;
            s->l1 = 0;
            p->n++;
            if ( p->n >= chunk_size ) break;
        }
        while(1);
    }
    else {
        do {
            if ( kseq_read(k1) < 0 ) break;
            if ( kseq_read(k2) < 0 ) break;
            trim_read_tail(k1->name.s, k1->name.l);
            trim_read_tail(k2->name.s, k2->name.l);
            
            if ( check_name(k1->name.s, k2->name.s) ) error("Inconsistance paired read names. %s vs %s.", k1->name.s, k2->name.s);
            //if ( k1->seq.l != k2->seq.l ) error("Inconsistant PE read length, %s.", k1->name.s);
            
            struct bseq *s;
            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct bseq));
            }
            s = &p->s[p->n];
            memset(s, 0, sizeof(*s));
            s->n0 = strdup(k1->name.s);
            s->s0 = strdup(k1->seq.s);
            s->q0 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l0 = k1->seq.l;
            //s->n1 = strdup(k2->name.s);
            s->s1 = strdup(k2->seq.s);
            s->q1 = k2->qual.l? strdup(k2->qual.s) : 0;
            s->l1 = k2->seq.l;
            
            p->n++;
            if ( p->n >= chunk_size ) break;            
        }
        while(1);
    }
    if ( p->n == 0 ) {
        bseq_pool_destroy(p);
        return NULL;
    }
return p;
}

struct fastq_handler *fastq_handler_init(const char *r1, const char *r2, int smart, int chunk_size)
{
    struct fastq_handler *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->r1 = strcmp(r1, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(r1, "r");
    if (h->r1 == NULL) error("Failed to open %s : %s.", r1, strerror(errno));
    h->k1 = kseq_init(h->r1);
    if (h->k1 == NULL) error("Failed to init stream. %s", r1);
    
    if (r2) {
        h->r2 = gzopen(r2, "r");    
        if (h->r2 == NULL) error("Failed to open %s: %s.", r2, strerror(errno));
        h->k2 = kseq_init(h->r2);
    }
    h->smart_pair = smart;
    h->chunk_size = chunk_size;
    return h;
}
void fastq_handler_destory(struct fastq_handler *h)
{
    kseq_destroy(h->k1);
    gzclose(h->r1);
    if ( h->k2 ) {
        kseq_destroy(h->k2);
        gzclose(h->r2);
    }
    free(h);
}
int fastq_handler_state(struct fastq_handler *h)
{
    if ( h == NULL ) return FH_NOT_ALLOC;
    if ( h->k1 == NULL ) return FH_NOT_INIT;
    if ( h->smart_pair ) return FH_SMART_PAIR;
    if ( h->k2 == NULL ) return FH_SE;
    return FH_PE;    
}

void *bseq_read(void *_h, void *opts)
{
    struct fastq_handler *h = (struct fastq_handler*)_h;

    int state = fastq_handler_state(h);
    struct bseq_pool *b = NULL;

    switch (state) {
        
        case FH_SE:
            b = bseq_read_core((kseq_t*)h->k1, NULL, h->chunk_size, 0);
            break;
            
        case FH_PE:
            b = bseq_read_core((kseq_t*)h->k1, (kseq_t*)h->k2, h->chunk_size, 1);
            break;
            
        case FH_SMART_PAIR:
            b = bseq_read_smart((kseq_t*)h->k1, h->chunk_size);
            break;
            
        case FH_NOT_ALLOC:
            error("The fastq handler is NOT allocated.");
            break;
        case FH_NOT_INIT:
            error("The fastq handler is NOT inited.");
            break;
        default:
            error("Unknown state, %d.", state);
            break;
    }
    if (b) b->opts = opts;
    return b;
}

int fastq_process(struct bseq_pool *pool, void *(*func)(void *data))
{
    int i;
    for ( i = 0; i < pool->n; ++i ) {
        struct bseq *s = &pool->s[i];
        s = func(s);
    }
    return 0;
}

#include "htslib/khash.h"

struct hval {
    int idx;
    int qual;
};
struct hvals {
    int n, m;
    struct hval *v;
};
KHASH_MAP_INIT_STR(key, struct hvals*)
//
// Process fastq block...
// stash
void bseq_pool_push(struct bseq *b, struct bseq_pool *p)
{
    assert(p);
    if (p->n == p->m) {
        p->m = p->m == 0 ? 10 : p->m<<1;
        p->s = realloc(p->s, sizeof(struct bseq)*p->m);               
    }
    struct bseq *c = &p->s[p->n++];
    memcpy(c, b, sizeof(struct bseq));
    //debug_print("c %d, b %d", c->l0, b->l0);
    free(b);
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

static char *reverse_seq(char *s, int l)
{
    char *r = malloc(sizeof(char)*l);
    int i;
    for (i = 0; i < l; ++i ) {
        switch(s[i]) {
            case 'A':
                r[l-i-1]='T'; break;
            case 'C':
                r[l-i-1]='G'; break;
            case 'G':
                r[l-i-1]='C'; break;
            case 'T':
                r[l-i-1]='A'; break;
            default:
                error("Unknown bases, %s",s);
        }
    }
    return r;
}
static int check_dup(struct bseq *r, struct bseq *q, int strand)
{
    if (strand == -1) error("Unknown strand.");

    kstring_t str = {0,0,0};
    kstring_t str1 = {0,0,0};
    if (r->l1 > 0 && q->l1 > 0) {
        if (strand == 1) {
            int l = q->l1 > r->l0 ? r->l0 : q->l1;
            kputsn(r->s0, l, &str);
            char *rs = reverse_seq(q->s1, q->l1);
            kputsn(rs, l, &str1);
            free(rs);

            // now to read 2

            l = q->l0 > r->l1 ? r->l1 : q->l0;
            kputsn(r->s1+r->l1-l, l, &str);
            char *rs1 = reverse_seq(q->s0, q->l0);
            kputsn(rs1+q->l1-1, l, &str1);
            free(rs1);
        }
        else {
            int l = q->l0 > r->l0 ? r->l0 : q->l0;
            kputsn(q->s0, l, &str);
            kputsn(q->s1, l, &str1);
            l = q->l1 > r->l1 ? r->l1 : q->l1;
            kputsn(r->s1+r->l1-l, l, &str);
            kputsn(q->s1+q->l1-l, l, &str1);
        }
    }
    else {
        if (strand == 1 && q->l0 != r->l0) return 0;
        int l = q->l0 < r->l0 ? q->l0 : r->l0;
        kputsn(r->s0, l, &str);
        kputsn(q->s0, l, &str1);
    }

    kputs("", &str);
    kputs("", &str1);

    assert(str.l == str1.l);
    if (levenshtein(str.s, str1.s, str.l) > 2) {
        free(str.s);
        free(str1.s);
        return 0;
    }

    free(str.s);
    free(str1.s);
    return 1; // on dup
}
#define DEDUP_SEED 16
// require all sequence should greater than 16bp, and equal length
// the function will destroy all the flags marked before
int bseq_pool_dedup(struct bseq_pool *p)
{
    kh_key_t *hash = kh_init(key);
    int n =0, m = 0;
    char **key = NULL;

    int i;
    int j;
    khint_t k;

    kstring_t seed={0,0,0};
    kstring_t rseed = {0,0,0};
    
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        if (b->l0 < DEDUP_SEED) error("Read length is too short. %d", b->l0);
            
        // todo:: convert mark to bits
        b->flag = 0; // reset all the mark
        int qual = 0;
        
        if (b->q0) 
            for (j = 0; j < b->l0; ++j) qual += b->q0[j]-33;
        int strand = -1;
        seed.l = 0; rseed.l = 0; // reset seed string
        kputsn(b->s0, DEDUP_SEED, &seed);
        kputs("", &seed);
        char *rs = b->l1 >0 ? b->s1 : b->s0;
        int l_rs = b->l1 >0 ? b->l1 : b->l0;
        for (j = 0; j < DEDUP_SEED; ++j) {
            switch(rs[l_rs-1-j]) {
                case 'A':
                    kputc('T', &rseed);
                    break;
                case 'C':
                    kputc('G', &rseed);
                    break;
                case 'G':
                    kputc('C', &rseed);
                    break;
                case 'T':
                    kputc('A', &rseed);
                    break;
                default:
                    error("Try to reverse unknown base, %s",rs);
            }
            kputs("", &rseed);            
        }
        
        k = kh_get(key, hash, seed.s);
        if (k == kh_end(hash)) {
            // check reverse then
            k = kh_get(key, hash, rseed.s);
            if (k == kh_end(hash)) goto push_to_index; // reverse do not also match
            strand = 1;
            goto check_dup_records;
            
        }
        else {
            strand = 0;
            goto check_dup_records;
        }
        if (0) {
          check_dup_records:
            do {
                struct hvals *v;
                v = kh_val(hash, k);
                for (j = 0; j < v->n; ++j) {
                    struct hval *v1 = &v->v[j];
                    struct bseq *r = &p->s[v1->idx];
                    if (check_dup(r, b, strand)) {
                        //debug_print("%s\t%s\t%d", r->s0, b->s0, v1->idx);
                        if (qual > v1->qual) {
                            // update
                            v1->idx = i;
                            v1->qual = qual;
                            r->flag = FQ_DUP; // set last as dup
                        }
                        else {
                            b->flag = FQ_DUP;
                        }
                        break;
                    }            
                }
                // no found, push to index      
                if (j == v->n) goto push_to_index;
            } while(0);
        }

        if (0) {
          push_to_index:
            // only push forward string to index
            k = kh_get(key, hash, seed.s);
            if (k == kh_end(hash)) {
                if (m == n) {
                    m = m== 0? 10 : m*2;
                    key = realloc(key, m*sizeof(char*));
                }
                key[n] = strdup(seed.s);
                int ret;
                k = kh_put(key, hash, key[n], &ret);
                struct hvals *v = malloc(sizeof(*v));
                memset(v, 0, sizeof(*v));
                v->m = 2; 
                v->v = malloc(sizeof(struct hval)*v->m);
                struct hval *v1 = &v->v[v->n++];
                v1->idx = i;
                v1->qual = qual;
                kh_val(hash, k) = v;                    
                n++; // increase key index
            }
            else {
                struct hvals *v;
                v = kh_val(hash, k);
                if (v->n == v->m) {
                    v->m = v->m*2; // v->m never equal 0
                    v->v = realloc(v->v, v->m*sizeof(struct hval));
                }
                struct hval *v1 = &v->v[v->n++];
                v1->idx = i;
                v1->qual = qual;
            }
        }
    }


    for (i = 0; i < n; ++i) {
        k = kh_get(key, hash, key[i]);
        if (k == kh_end(hash)) error("Not in the key %s", key[i]);
        struct hvals *v = kh_val(hash, k);
        free(v->v);
        free(v);
        free(key[i]);
    }
    free(key);
    kh_destroy(key, hash);
    return 0;   
}
