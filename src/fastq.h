#ifndef FASTQ_H
#define FASTQ_H

#include "utils.h"
#include<zlib.h>

struct qc_report {
    uint64_t all_fragments;
    uint64_t qc_failed;
    uint64_t unknown_barcodes;
};

#define FQ_PASS     0
#define FQ_QC_FAIL  1
#define FQ_BC_FAIL  2

struct bseq {
    int flag; // flag for skip reasons
    // read 1
    char *n0;
    char *s0, *q0;
    int l0;
    // read 2
    char *s1, *q1;
    int l1;

    void *data; // extend data, should be freed manually
};

struct bseq_pool {
    int n, m;
    struct bseq *s;
    void *opts; // used to point thread safe structure
};

struct fastq_handler {
    gzFile r1;
    gzFile r2;
    void *k1;
    void *k2;
    int smart_pair;
    int chunk_size;
};

#define FH_SE 1
#define FH_PE 2
#define FH_SMART_PAIR 3
#define FH_NOT_ALLOC 4
#define FH_NOT_INIT 5

extern int check_name(char *s1, char *s2);

struct bseq_pool *bseq_pool_init();

void bseq_pool_clean(struct bseq_pool *p);

void bseq_pool_destroy(struct bseq_pool *p);

// fastq handler must be inited before call bseq_read
void *bseq_read(void *h, void *opts);

extern struct fastq_handler *fastq_handler_init(const char *r1, const char *r2, int smart, int chunk_size);

extern int fastq_handler_state(struct fastq_handler*);

extern void fastq_handler_destory(struct fastq_handler *h);


#endif
