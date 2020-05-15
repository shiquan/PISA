#ifndef FASTQ_H
#define FASTQ_H

#include "utils.h"
#include "dict.h"
#include <zlib.h>

struct bseq_core {
    int length;
    char *seq;
    char *qual;
}
    
struct bseq {
    int flag; 
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


#endif
