#ifndef BED_H
#define BED_H

#include "utils.h"
#include "region_index.h"
#include "dict.h"

struct bed {    
    int seqname;
    int name;
    int start;
    int end;
    // float score; // if column 5 exists
    int strand;  // 0 on forward, 1 on backword    
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

struct bed_spec *bed_spec_init();
void bed_spec_destroy(struct bed_spec *B);
struct bed_spec *bed_read(const char *fname);
struct region_itr *bed_query(struct bed_spec *B, char *name, int start, int end);

#endif
