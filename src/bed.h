#ifndef BED_H
#define BED_H

#include "utils.h"
#include "region_index.h"
#include "dict.h"
#include "htslib/kstring.h"

#define BED_STRAND_FWD 0
#define BED_STRAND_REV 1
#define BED_STRAND_UNK -1

struct bed {    
    int seqname;
    int name;
    int start;
    int end;
    // float score; // if column 5 exists
    int strand;  // 0 on forward, 1 on backword
    void *data;
};

struct _ctg_idx;

struct bed_spec {
    struct dict *seqname;
    struct dict *name;    

    struct bed_idx *idx;
    struct _ctg_idx *ctg;
    int n,m;
    struct bed *bed;
};

struct var {
    kstring_t *ref;
    kstring_t *alt;
};

struct bed_spec *bed_spec_init();
void bed_spec_destroy(struct bed_spec *B);
struct bed_spec *bed_read(const char *fname);
struct region_itr *bed_query(const struct bed_spec *B, char *name, int start, int end);
int bed_check_overlap(const struct bed_spec *B, char *name, int start, int end);


struct bed_spec *bed_read_vcf(const char *fn);
void bed_spec_var_destroy(struct bed_spec *B);
#endif
