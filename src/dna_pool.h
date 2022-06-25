#ifndef DNA_POOL_H
#define DNA_POOL_H
#include "utils.h"

struct PISA_dna {
    union {
        int idx;
        uint8_t *dna;
    };
    union {
        int count;
        void *data;
    };
};

struct PISA_dna_pool {
    struct PISA_dna *data;
    int l, m;
    int len; // dna length, all DNAs in the pool should be equal length
};

void PISA_dna_destroy(struct PISA_dna_pool *p);
void PISA_idx_destroy(struct PISA_dna_pool *p);

struct PISA_dna_pool *PISA_dna_pool_init();
//struct PISA_dna_pool *PISA_idx_init();

//int PISA_dna_query(struct PISA_dna_pool *p, const char *seq);
//int PISA_idx_query(struct PISA_dna_pool *p, const int idx);

struct PISA_dna *PISA_dna_push(struct PISA_dna_pool *p, const char *seq);
struct PISA_dna *PISA_dna_push1(struct PISA_dna_pool *p, const char *seq, void *data);
struct PISA_dna *PISA_dna_query(struct PISA_dna_pool *p, const char *seq);

struct PISA_dna *PISA_idx_push(struct PISA_dna_pool *p, const int idx);
struct PISA_dna *PISA_idx_push1(struct PISA_dna_pool *p, const int idx, void *data);
struct PISA_dna *PISA_idx_query(struct PISA_dna_pool *p, const int idx);

struct PISA_dna_pool *PISA_pool_merge(struct PISA_dna_pool *p1,struct PISA_dna_pool *p2);
                                      

#endif
