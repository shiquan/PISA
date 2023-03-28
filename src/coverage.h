#ifndef COVERAGE_H
#define COVERAGE_H

#include "utils.h"
#include "dict.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"

struct depth {
    int pos;
    
    int dep1; // unknown strand or forward
    int dep2; // reverse strand
    struct dict *bc1;
    struct dict *bc2;

    // int depth;
    // int strand;
    int id;
    // struct dict *bc;
    struct depth *next;
    struct depth *before;

    // a subtree for different ids
    struct depth *left;
    struct depth *right;
};

struct depth *depth_init();

void depth_destroy(struct depth *d);

struct depth* bam2depth(const hts_idx_t *idx, const int tid, const int start, const int end,
                        const int strand,
                        htsFile *fp,
                        const int mapq_thres,
                        const int ignore_strand,
                        struct dict *bc,
                        const char *bc_tag,
                        const char *umi_tag,
                        const int split_by_tag,
                        const int alias_tag,
                        const int *alias_idx,
                        int fix_barcodes
    );

#endif
