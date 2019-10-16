#ifndef FASTQ_TD
#define FASTQ_TD
#include <stdio.h>
#include <stdlib.h>
#include "dict.h"

struct read {
    char *name;
    int l0, l1;
    char *s0;
    char *q0;
    char *s1;
    char *q1;
};

struct read_block {
    char *name;
    struct read *b;
    int n, m;
    int pair_mode;
};

struct thread_dat {
    int n, m;
    struct read_block *rb;
};

void read_block_push(struct read_block *b, struct read *r);
struct read_block *read_block_copy(struct read_block *r);
void read_block_clear(struct read_block *b);

extern void thread_dat_destroy(struct thread_dat *td);
extern struct thread_dat *read_thread_dat(FILE *fp, struct dict *dict);
    
#endif
