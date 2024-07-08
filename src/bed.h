#ifndef BED_H
#define BED_H

#include "utils.h"
#include "region_index.h"
#include "dict.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "gtf.h"

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

void bed_copy(struct bed *a, struct bed *b);

#define BAT_COUNT    22

// definition of BED annotation types
#define BAT_UNKNOWN          0
// promoter come before UTR considering region behind TSS may also play as promoter,
// but maybe this is not a perfect design, may change in the future. SQ 19/04/2024
#define BAT_PROMOTER         1
#define BAT_UTR3             2 
#define BAT_UTR5             3 
#define BAT_EXON             4 
#define BAT_MULTIEXONS       5 
#define BAT_EXONINTRON       6
#define BAT_WHOLEGENE        7 
#define BAT_MULTIGENES       8 
#define BAT_INTRON           9 
#define BAT_ANTISENSEUTR3    10
#define BAT_ANTISENSEUTR5    11
#define BAT_ANTISENSEEXON    12
#define BAT_ANTISENSEINTRON  13
#define BAT_ANTISENSECOMPLEX 14
#define BAT_FLANK            15
#define BAT_UPSTREAM         16
#define BAT_DOWNSTREAM       17
#define BAT_ANTISENSEUP      18
#define BAT_ANTISENSEDOWN    19
#define BAT_INTERGENIC       20
#define BAT_UNKNOWNCHRS      21

extern const char *bed_typename(int type);

struct bed_ext {
    int n;
    char **genes;

    int type;

    // distance to nearby gene, 0 for enclosed region
    // int distance; 
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
    void *ext;
};

struct bed_spec *bed_spec_init();

void bed_spec_ext_destroy(struct bed_spec *B);
void bed_spec_destroy(struct bed_spec *B);

struct bed_spec *bed_spec_dup(struct bed_spec *B0);
struct bed_spec *bed_read0(struct bed_spec *B, const char *fname);
struct bed_spec *bed_read(const char *fname);
void bed_build_index(struct bed_spec *B);
void bed_spec_sort(struct bed_spec *B);
// start is 0 based
struct region_itr *bed_query(const struct bed_spec *B, char *name, int start, int end, int strand);
int bed_check_overlap(const struct bed_spec *B, char *name, int start, int end, int strand);
char* bed_seqname(struct bed_spec *B, int id);
int bed_name2id(struct bed_spec *B, char *name);
int bed_spec_push0(struct bed_spec *B, const char *seqname, int start, int end, int strand, const char *name, void *ext);
int bed_spec_push1(struct bed_spec *B, int seqname, int start, int end, int strand, int name, void *ext);
int bed_spec_push(struct bed_spec *B, struct bed *bed);
struct bed_spec *bed_read_vcf(const char *fn);

void bed_spec_merge0(struct bed_spec *B, int strand, int check_name);
void bed_spec_merge1(struct bed_spec *B, int strand, int up, int down, int min_length, int check_name);
void bed_spec_merge2(struct bed_spec *B, int strand, int gap, int min_length, int check_name);

void bed_spec_var_destroy(struct bed_spec *B);
void bed_spec_write0(struct bed_spec *B, FILE *out, int ext, int gene_as_name);
void bed_spec_write(struct bed_spec *B, const char *fn, int ext, int gene_as_name);
void bed_spec_seqname_from_bam(struct bed_spec *B, bam_hdr_t *hdr);
void bed_spec_dedup(struct bed_spec *B, int check_name);
struct bed_spec *bed_spec_flatten(struct bed_spec *B, int offset);
struct bed_spec *gtf2bed(struct gtf_spec *G, struct region_itr *itr, int level, int name_level, int offset);
// bed anno
struct anno0 {
    struct gtf *g;
    int type;
};

struct anno0 *anno_bed_core(const char *name, int start, int end, int strand, struct gtf_spec *G, int *n, int down, int up);

void anno_bed_cleanup();

#endif
