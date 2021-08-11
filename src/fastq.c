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
    for(n = 0; n < l1; ++n, ++s1, ++s2)
        if (*s1 != *s2) return 1;
    
    return 0;
}

struct bseq_pool *bseq_pool_init()
{
    struct bseq_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    return p;
}
void bseq_clean(struct bseq *b)
{
    if (b->n0.m) free(b->n0.s);
    if (b->s0.m) free(b->s0.s);
    if (b->q0.m) free(b->q0.s);
    if (b->s1.m) free(b->s1.s);
    if (b->q1.m) free(b->q1.s);
}
void bseq_destroy(struct bseq *b)
{
    bseq_clean(b);
    free(b);
}
void bseq_pool_clean(struct bseq_pool *p)
{
    int i;
    for ( i = 0; i < p->n; ++i ) {
        struct bseq *b = &p->s[i];
        bseq_clean(b);
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

kstring_t *kstr_init()
{
    kstring_t *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(*s));
    return s;
}
void bseq_init(struct bseq *b)
{
    memset(b, 0, sizeof(*b));
}

int kstr_copy(kstring_t *a, kstring_t *b)
{
    if (b->l == 0 || b->s == 0) {
        warnings("Try to copy an empty string.");
        return 1;
    }
    a->l = 0;
    kputsn(b->s, b->l, a);
    kputs("",a);

    return 0;
}

static struct bseq_pool *fastq_read_smart(struct fastq_handler *h, int chunk_size)
{
    struct bseq_pool *p = bseq_pool_init();
    int size = 0;
    int ret1= -1;
    do {
        
        ret1 = kseq_read(h->k1);
    
        if (ret1 < 0) { // come to the end of file
            if (h->n_file > 1 && h->curr < h->n_file) {
                gzclose(h->r1);
                kseq_destroy(h->k1);
                h->r1 = gzopen(h->read_1[h->curr], "r");
                h->k1 = kseq_init(h->r1);
                if (kseq_read(h->k1) < 0) error("Empty record ? %s", h->read_1[h->curr]);
            }
            else break;
            h->curr++;
        }

        kseq_t *ks = h->k1;
        
        struct bseq *s;
        if (p->n >= p->m ) {
            p->m = p->m ? p->m*2 : 256;
            p->s = realloc(p->s, p->m*sizeof(struct bseq));
        }
        s = &p->s[p->n];

        trim_read_tail(ks->name.s, ks->name.l);

        bseq_init(s);
        
        kstr_copy(&s->n0, &ks->name);
        kstr_copy(&s->s0, &ks->seq);
        kstr_copy(&s->q0, &ks->qual);

        if ( kseq_read(ks) < 0 ) error("Truncated input.");

        trim_read_tail(ks->name.s, ks->name.l);
        
        if ( check_name(s->n0.s, ks->name.s) ) error("Inconsistance paired read names. %s vs %s.", s->n0.s, ks->name.s);

        kstr_copy(&s->s1, &ks->seq);
        kstr_copy(&s->q1, &ks->qual);
        
        size += s->s0.l;
        size += s->s1.l;
        p->n++;
        if ( size >= chunk_size ) break;
    } while (1);
    
    if ( p->n == 0 ) {
        bseq_pool_destroy(p);
        return NULL;
   }
    return p;
}
static struct bseq_pool *fastq_read_core(struct fastq_handler *h, int chunk_size, int pe)
{
    // k1 and k2 already load one record when come here
    struct bseq_pool *p = bseq_pool_init();
    int ret1, ret2 = -1;
    
    if ( pe == 0 ) {
        do {
            ret1 = kseq_read(h->k1);
            
            if (ret1 < 0) { // come to the end of file
                if (h->n_file > 1 && h->curr < h->n_file) {
                    gzclose(h->r1);
                    kseq_destroy(h->k1);
                    h->r1 = gzopen(h->read_1[h->curr], "r");
                    h->k1 = kseq_init(h->r1);
                    if (kseq_read(h->k1) < 0) error("Empty record ? %s", h->read_1[h->curr]);
                }
                else break;
                h->curr++;
            }

            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct bseq));
            }
            struct bseq *s = &p->s[p->n];
            kseq_t *k1 = h->k1;
            trim_read_tail(k1->name.s, k1->name.l);
            bseq_init(s);
            kstr_copy(&s->n0, &k1->name);
            kstr_copy(&s->s0, &k1->seq);
            kstr_copy(&s->q0, &k1->qual);
            p->n++;
            if ( p->n >= chunk_size ) break;
        }
        while(1);
    }
    else {
        do {            
            ret1 = kseq_read(h->k1);
            ret2 = kseq_read(h->k2);
            if (ret1 < 0) { // come to the end of file
                if (ret2 >=0) error("Inconsistant input fastq records.");
                if (h->n_file > 1 && h->curr < h->n_file) {
                    gzclose(h->r1);
                    kseq_destroy(h->k1);
                    h->r1 = gzopen(h->read_1[h->curr], "r");
                    h->k1 = kseq_init(h->r1);
                    if (kseq_read(h->k1) < 0) error("Empty record ? %s", h->read_1[h->curr]);
                    if (h->r2) {
                        gzclose(h->r2);
                        kseq_destroy(h->k2);
                        h->r2 = gzopen(h->read_2[h->curr], "r");
                        h->k2 = kseq_init(h->r2);
                        if (kseq_read(h->k2) < 0) error("Empty record ? %s", h->read_2[h->curr]);
                    }
                    h->curr++;
                }
                else break;
            }
            kseq_t *k1 = h->k1;
            kseq_t *k2 = h->k2;
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
            bseq_init(s);
            kstr_copy(&s->n0, &k1->name);
            kstr_copy(&s->s0, &k1->seq);
            kstr_copy(&s->q0, &k1->qual);
            kstr_copy(&s->s1, &k2->seq);
            kstr_copy(&s->q1, &k2->qual);
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
static char **split_multi_files(const char *fname, int *n)
{
    kstring_t str = {0,0,0};
    kputs(fname, &str);
    int *s = ksplit(&str, ',', n);
    if (*n == 1 || s == 0) {
        free(str.s);
        if (s) free(s);
        return NULL;
    }
    
    char **paths = malloc(*n*sizeof(char*));
    int i;
    for (i = 0; i < *n; ++i) paths[i] = strdup(str.s+s[i]);
    free(str.s);
    free(s);
    return paths;
}

struct fastq_handler *fastq_handler_init(const char *r1, const char *r2, int smart, int chunk_size)
{
    struct fastq_handler *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    int n1, n2;
    h->read_1 = split_multi_files(r1, &n1);
    if (r2 != NULL) {
        h->read_2 = split_multi_files(r2, &n2);
        if (n1 != n2) error("Unpaired input fastqs.");
    }
    assert(n1 > 0);
    h->n_file = n1;
    h->curr = 1;
    h->smart_pair = smart;
    h->chunk_size = chunk_size;

    if (n1 == 1) {
        h->r1 = strcmp(r1, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(r1, "r");
        if (h->r1 == NULL) error("Failed to open %s : %s.", r1, strerror(errno));
        h->k1 = kseq_init(h->r1);
        if (h->k1 == NULL) error("Failed to init stream. %s", r1);
        
        if (r2) {
            h->r2 = gzopen(r2, "r");    
            if (h->r2 == NULL) error("Failed to open %s: %s.", r2, strerror(errno));
            h->k2 = kseq_init(h->r2);
        }
    }
    else {
        h->r1 = gzopen(h->read_1[0], "r");
        h->k1 = kseq_init(h->r1);
        if (r2) {
            h->r2 = gzopen(h->read_2[0], "r");
            h->k2 = kseq_init(h->r2);
        }
    }
    
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
    if (h->n_file > 1) {
        int i;
        for (i = 0; i < h->n_file;++i) {
            free(h->read_1[i]);
            if (h->read_2) free(h->read_2[i]);
        }
        free(h->read_1);
        if (h->read_2) free(h->read_2);
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
void *fastq_read(void *_h, void *opts)
{
    struct fastq_handler *h = (struct fastq_handler*)_h;

    int state = fastq_handler_state(h);

    struct bseq_pool *b;
    
    switch(state) {
        case FH_SE:
            b = fastq_read_core(h, h->chunk_size, 0);
            break;
            
        case FH_PE:
            b = fastq_read_core(h, h->chunk_size, 1);
            break;
            
        case FH_SMART_PAIR:
            b = fastq_read_smart(h, h->chunk_size);
            break;

        case FH_NOT_ALLOC:
            error("The fastq handler is NOT allocated.");
            break;

        case FH_NOT_INIT:
            error("The fastq handler is NOT inited.");
            break;

        default:
            error("Unknown state");
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
    // free(b);
}

size_t hamming_n(const char *a, const size_t length, const char *b, const size_t bLength) {
    if (a == b) return 0;
    if (length == 0) return bLength;
    if (bLength == 0) return length;
    int i;
    int dist = 0;
    for (i = 0; i < length && i < bLength; ++i)
        if (a[i] != b[i]) dist++;
    
    if (i < length) dist += length-i;
    if (i < bLength) dist += bLength-i;
    return dist;
}
// credit to https://github.com/wooorm/levenshtein.c
size_t levenshtein_n(const char *a, const size_t length, const char *b, const size_t bLength) {
  // Shortcut optimizations / degenerate cases.
  if (a == b) {
    return 0;
  }

  if (length == 0) {
    return bLength;
  }

  if (bLength == 0) {
    return length;
  }

  size_t *cache = calloc(length, sizeof(size_t));
  size_t index = 0;
  size_t bIndex = 0;
  size_t distance;
  size_t bDistance;
  size_t result;
  char code;

  // initialize the vector.
  while (index < length) {
    cache[index] = index + 1;
    index++;
  }

  // Loop.
  while (bIndex < bLength) {
    code = b[bIndex];
    result = distance = bIndex++;
    index = SIZE_MAX;

    while (++index < length) {
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
            case 'N':
                r[l-i-1]='N'; break;
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
    if (r->s0.l > 0 && q->s1.l > 0) {
        if (strand == 1) {
            int l = q->s1.l > r->s0.l ? r->s0.l : q->s1.l;
            kputsn(r->s0.s, l, &str);
            char *rs = reverse_seq(q->s1.s, q->s1.l);
            kputsn(rs, l, &str1);
            free(rs);

            // now to read 2

            l = q->s0.l > r->s1.l ? r->s1.l : q->s0.l;
            kputsn(r->s1.s+r->s1.l-l, l, &str);
            char *rs1 = reverse_seq(q->s0.s, q->s0.l);
            kputsn(rs1+q->s1.l-1, l, &str1);
            free(rs1);
        }
        else {
            int l = q->s0.l > r->s0.l ? r->s0.l : q->s0.l;
            kputsn(q->s0.s, l, &str);
            kputsn(q->s1.s, l, &str1);
            l = q->s1.l > r->s1.l ? r->s1.l : q->s1.l;
            kputsn(r->s1.s+r->s1.l-l, l, &str);
            kputsn(q->s1.s+q->s1.l-l, l, &str1);
        }
    }
    else {
        if (strand == 1 && q->s0.l != r->s0.l) return 0;
        int l = q->s0.l < r->s0.l ? q->s0.l : r->s0.l;
        kputsn(r->s0.s, l, &str);
        kputsn(q->s0.s, l, &str1);
    }

    kputs("", &str);
    kputs("", &str1);

    assert(str.l == str1.l);
    //int score = levenshtein_n(str.s, str.l, str1.s, str1.l);
    int score = hamming_n(str.s, str.l, str1.s, str1.l);
    if (score > 2) {
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
    if (p->n == 1) return 0;
    
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
        if (b->s0.l < DEDUP_SEED) error("Read length is too short. %zu", b->s0.l);
            
        // todo:: convert mark to bits
        b->flag = 0; // reset all the mark
        int qual = 0;
        
        if (b->q0.l) 
            for (j = 0; j < b->s0.l; ++j) qual += b->q0.s[j]-33;
        int strand = -1;
        seed.l = 0; rseed.l = 0; // reset seed string
        kputsn(b->s0.s, DEDUP_SEED, &seed);
        kputs("", &seed);
        char *rs = b->s1.l >0 ? b->s1.s : b->s0.s;
        int l_rs = b->s1.l >0 ? b->s1.l : b->s0.l;
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
                case 'N':
                    kputc('N', &rseed);
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

    // debug_print("Start debug");
    for (i = 0; i < n; ++i) {
        k = kh_get(key, hash, key[i]);
        if (k == kh_end(hash)) error("Not in the key %s", key[i]);
        struct hvals *v = kh_val(hash, k);
        // debug_print("%s", key[i]);
        free(v->v);
        free(v);
        free(key[i]);
        kh_del(key,hash, k);
    }
    free(seed.s);
    free(rseed.s);
    free(key);
    kh_destroy(key, hash);
    return 0;   
}
/* TODO: improve buffer performance
#define MIN_BUFFER      1     // 1M
#define DEFAULT_BUFFER  1000  // 1G
#define MAX_BUFFER      10000 // 10G

struct fastq_buffer {
    uint32_t n, m;
    struct bseq *s;
    int paired;
    uint32_t l_buf, m_buf;
    uint8_t *buf;
};

struct fastq_handler {
    gzFile read_1;
    gzFile read_2;
    uint8_t *buf;
}
*/
