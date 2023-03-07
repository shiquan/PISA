#ifndef RD_ANNO_HEADER
#define RD_ANNO_HEADER
#include "utils.h"

// todo: merge with bed annotation types
enum exon_type {
    type_unknown = 0,  // unknown type, init state
    type_exon,     // read full covered in exon
    type_intron,        // read full covered in intron, with same strand of gene
    type_exon_intron,   // read cover exon and nearby intron
    type_antisense,     // read map on antisense
    type_splice,        // junction read map two or more exome
    type_ambiguous,     // junction read map to isoform(s) but skip some isoforms between, or map to intron
    type_intergenic,
    type_antisense_intron,
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

enum exon_type RE_type_map(uint8_t c);
const char *RE_tag_name(int i);
const char *exon_type_name(int i);
#endif
