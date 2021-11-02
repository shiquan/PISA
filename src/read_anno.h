#ifndef RD_ANNO_HEADER
#define RD_ANNO_HEADER
#include "utils.h"

enum exon_type {
    type_unknown = 0,  // unknown type, init state
    type_exon,     // read full covered in exon
    type_intron,        // read full covered in intron, with same strand of gene
    type_exon_intron,   // read cover exon and nearby intron
    type_antisense,     // read map on antisense
    type_splice,        // junction read map two or more exome
    type_ambiguous,     // junction read map to isoform(s) but skip some isoforms between, or map to intron
    type_intergenic,
};

static char *exon_type_names[] = {
    "Unknown", "Exon", "Intron", "ExonIntron", "Antisense", "Splice", "Ambiguous", "Intergenic"
};


static char *RE_tags[] = {
    "U", "E", "N", "C", "A", "S", "V", "I",
};

struct trans_type {
    int trans_id;
    enum exon_type type;
};

struct gene_type {
    int gene_id;
    int gene_name;
    enum exon_type type; // main type
    int n, m;
    struct trans_type *a;
};

struct gtf_anno_type {
    enum exon_type type;
    int n, m;
    struct gene_type *a;
};

static enum exon_type RE_map[256] = {
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  4,  0,  3,  0,  1,  0,  0,  0,  7,  0,  0,  0,  0,  2,  0, 
     0,  0,  0,  5,  0,  0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
};

static enum exon_type RE_type_map(char c) {
    return RE_map[(uint8_t)c];
}

#endif
