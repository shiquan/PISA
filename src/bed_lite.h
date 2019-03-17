#ifndef BED_HEADER
#define BED_HEADER

#include "utils.h"

struct bed {
    int start;
    int end;
    char *name;
    void *data;
};

struct bed_chr {
    int id;
    int n, m;
    struct bed *b;
};

struct bedaux {
    int n, m;
    char **names;
    struct bed_chr *c;
    void *name2id;
};

extern struct bedaux *bed_read(const char *fn);
extern int bed_select_chrom(struct bedaux *bed, const char *chr);
extern void bed_destroy(struct bedaux*);

#endif
