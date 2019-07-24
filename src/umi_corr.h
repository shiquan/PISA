#ifndef UMI_CORR_H
#define UMI_CORR_H

#include <stdlib.h>

struct block_umi;

struct corr_tag {
    int n, m;
    struct block_umi *umi;
    void *dict;
};

struct corr_tag *corr_tag_build();

void corr_tag_destory(struct corr_tag *C);

void corr_tag_push(struct corr_tag *C, char const *n, char const *u);

char *corr_tag_retrieve(struct corr_tag *C, char const *n, char const *u);

void corr_tag(struct corr_tag *C);

#endif
