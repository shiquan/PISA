#ifndef FASTQ_H
#define FASTQ_H

#include "utils.h"
#include "dict.h"
#include <zlib.h>
/*
struct qc_report {
    uint64_t all_fragments;
    uint64_t qc_failed;
    uint64_t unknown_barcodes;
};

#define FQ_PASS     0
#define FQ_QC_FAIL  1
#define FQ_BC_FAIL  2
#define FQ_DUP      3
*/
struct bseq_core {
    int length;
    char *seq;
    char *qual;
}
    
struct bseq {
    int flag; //reserved flag    
    int n;
    struct bseq_core *b;
    struct dict *extend_tags;
    uint8_t *data;
    int l_data;
    int m_data;
};

struct input;

struct bseq_pool {
    int n, m;
    struct bseq *s;
    void *data;
};

struct fastq_handler {
    int n_file;
    struct input *input;
    int input_smart_pair;
    struct bseq *buf;
};

extern struct fastq_handler *fastq_handler_init(char **input_fname, int n);
void fastq_handler_destroy(struct fastq_handler *fq);

struct bseq_pool *fastq_read(struct fastq_handler *fq, int n_record, int max_mem);

void bseq_pool_destroy(struct bseq_pool *p);
void bseq_pool_push(struct bseq_pool *p, struct bseq *b);
void bseq_destroy(struct bseq *b);
int bseq_pool_dedup(struct bseq_pool *p);

void *fastq_tag_value(struct bseq *b, const char *tag);
char *fastq_tags(struct bseq *b, struct dict *tags);

void fastq_tag_push(struct bseq *b, const char *tag, int type, void *data);

char *fastq_select_seq(struct bseq *b, int rd, int start, int end);
char *fastq_select_qual(struct bseq *b, int rd, int start, int end);
int fastq_mean_qual(struct bseq *b);

char *compact_long_DNA(char *seq);

/*
#define FH_SE 1
#define FH_PE 2
#define FH_SMART_PAIR 3
#define FH_NOT_ALLOC 4
#define FH_NOT_INIT 5

extern int check_name(char *s1, char *s2);

struct bseq_pool *bseq_pool_init();

void bseq_pool_clean(struct bseq_pool *p);

void bseq_pool_destroy(struct bseq_pool *p);

// fastq handler must be inited before call fastq_read
void *fastq_read(void *h, void *opts);

extern int fastq_handler_state(struct fastq_handler*);

extern void fastq_handler_destory(struct fastq_handler *h);
extern void bseq_pool_push(struct bseq *b, struct bseq_pool *p);
//extern int levenshtein(char *a, char *b, int l);
extern int bseq_pool_dedup(struct bseq_pool *p);

extern size_t levenshtein_n(const char *a, const size_t length, const char *b, const size_t bLength);
#endif
*/
