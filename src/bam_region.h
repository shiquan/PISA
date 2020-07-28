#ifndef BAM_REGION_H
#define BAM_REGION_H

#include "utils.h"
#include "bed.h"
#include "region_index.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

struct bam_region_itr {
    struct bed_spec *bed;
    hts_idx_t *idx;
    hts_itr_t *itr;
    int i_bed;
};

struct bam_region_itr *bri_init(htsFile *fp, const char *fn, struct bed_spec *bed);
int bri_read(htsFile *fp, bam_hdr_t *hdr, bam1_t *b, struct bam_region_itr *bri);
void bri_destroy(struct bam_region_itr *bri);
    

#endif
