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
