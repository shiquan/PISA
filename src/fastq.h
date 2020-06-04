#ifndef FASTQ_H
#define FASTQ_H

#include "utils.h"
#include "dict.h"
#include <zlib.h>

struct fastq0 {
    int length;
    char *seq;
    char *qual;
}
    
struct fastq {
    int flag; 
    int n_fastq;
    struct fastq0 *fastq;
    struct dict *extend_tags;
    int n_idx, m_idx;
    int *idx;
    uint8_t *data;
    int l_data;
    int m_data;
};

struct input;

struct fastq_pool {
    int n, m;
    struct fastq *fastq;
    void *data;
};

struct fastq_spec {
    int n_file;
    struct input *input;
    int input_smart_pair;
    struct fastq *buf; // temp buffer
};

extern struct fastq_spec *fastq_spec_init(char **input_fname, int n);
void fastq_init_destroy(struct fastq_init *fq);

struct fastq_pool *fastq_read(struct fastq_spec *fq, int n_record, int max_mem);

void fastq_pool_destroy(struct fastq_pool *p);
void fastq_pool_push(struct fastq_pool *p, struct fastq *b);
void fastq_destroy(struct fastq *b);
int fastq_pool_dedup(struct fastq_pool *p);

void *fastq_tag_value(struct fastq *b, const char *tag);
char *fastq_tags(struct fastq *b, struct dict *tags);

void fastq_tag_push(struct fastq *b, const char *tag, int type, void *data);

char *fastq_select_seq(struct fastq *b, int rd, int start, int end);
char *fastq_select_qual(struct fastq *b, int rd, int start, int end);
int fastq_mean_qual(struct fastq *b);

char *compact_long_DNA(char *seq);


#endif
