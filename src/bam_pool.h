#ifndef BAM_POOL_H
#define BAM_POOL_H

#include "htslib/hts.h"
#include "htslib/sam.h"

struct bam_pool {
    int n, m;
    bam1_t *bam;
};

extern struct bam_pool *bam_pool_create();
extern struct bam_pool *bam_pool_init(int size);
extern void bam_read_pool(struct bam_pool *p, htsFile *fp, bam_hdr_t *h, int chunk_size);
extern void bam_pool_destory(struct bam_pool *p);

#endif
