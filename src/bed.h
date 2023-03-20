#ifndef BED_H
#define BED_H

#include "utils.h"
#include "region_index.h"
#include "dict.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

#define BED_STRAND_FWD 0
#define BED_STRAND_REV 1
#define BED_STRAND_UNK -1
#define BED_STRAND_IGN -1

#define strand_is_minus(a) (a)==1

struct bed {    
    int seqname;
    int name;
    int start; // 0 base
    int end;   // 1 base
    // float score; // if column 5 exists
    int strand;  // 0 on forward, 1 on backword, -1 on unset or ignore
    void *data;
};

#define BAT_COUNT    15

// definition of BED annotation types
#define BAT_UNKNOWN          0
#define BAT_MULTIGENES       1
#define BAT_WHOLEGENE        2
#define BAT_UTR3             3
#define BAT_UTR5             4
#define BAT_EXON             5
#define BAT_MULTIEXONS       6
#define BAT_EXONINTRON       7
#define BAT_INTRON           8
#define BAT_ANTISENSEUTR3    9
#define BAT_ANTISENSEUTR5    10
#define BAT_ANTISENSEEXON    11
#define BAT_ANTISENSEINTRON  12
#define BAT_ANTISENSECOMPLEX 13
#define BAT_INTERGENIC       14

extern const char *bed_typename(int type);

struct bed_ext {
    int n;
    char **genes;

    int type;

    // distance to nearby gene, 0 for enclosed region
    int distance; 
};

extern struct bed_ext *bed_ext_init();

struct _ctg_idx;

struct bed_spec {
    struct dict *seqname;
    struct dict *name;    

    struct bed_idx *idx;
    struct _ctg_idx {
        int offset;
        int idx;
    } *ctg;
    // struct _ctg_idx *ctg;
    int n,m;
    struct bed *bed;
};

struct var {
    kstring_t *ref;
    kstring_t *alt;
};

struct bed_spec *bed_spec_init();

void bed_spec_ext_destroy(struct bed_spec *B);
void bed_spec_destroy(struct bed_spec *B);

struct bed_spec *bed_read0(struct bed_spec *B, const char *fname);
struct bed_spec *bed_read(const char *fname);

// start is 0 based
struct region_itr *bed_query(const struct bed_spec *B, char *name, int start, int end, int strand);
int bed_check_overlap(const struct bed_spec *B, char *name, int start, int end, int strand);
char* bed_seqname(struct bed_spec *B, int id);
int bed_name2id(struct bed_spec *B, char *name);
int bed_spec_push(struct bed_spec *B, struct bed *bed);
struct bed_spec *bed_read_vcf(const char *fn);

void bed_spec_merge0(struct bed_spec *B, int strand, int check_name);
void bed_spec_merge1(struct bed_spec *B, int strand, int up, int down, int min_length, int check_name);
void bed_spec_merge2(struct bed_spec *B, int strand, int gap, int min_length, int check_name);

void bed_spec_var_destroy(struct bed_spec *B);
void bed_spec_write0(struct bed_spec *B, FILE *out, int ext);
void bed_spec_write(struct bed_spec *B, const char *fn, int ext);
void bed_spec_seqname_from_bam(struct bed_spec *B, bam_hdr_t *hdr);
    
#endif
