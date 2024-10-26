#ifndef RD_ANNO_HEADER
#define RD_ANNO_HEADER
#include "utils.h"
#include "gtf.h"

// todo: merge with bed annotation types
enum exon_type {
    type_unknown = 0,  // unknown type, init state
    type_splice,        // junction read map two or more exome
    type_exon,     // read full covered in exon
    type_exon_intron,   // read cover exon and nearby intron
    type_intron_retain,
    type_intron,        // read full covered in intron, with same strand of gene
    type_exclude,       // exclude exons
    type_antisense,     // read map on antisense
    type_antisense_intron,
    type_ambiguous,     // junction read map to isoform(s) but skip some isoforms between, or map to intron
    type_intergenic,    
};

struct trans_type {
    int trans_id;
    enum exon_type type;
    int n_exon;
    int m_exon;
    struct gtf **exon;
    int n_exclude;
    int m_exclude;
    struct gtf **exl;
};

struct gene_type {
    int gene_id;
    int gene_name;
    enum exon_type type; // main type
    int n, m;
    struct trans_type *a;

    int n_flatten;
    int m_flatten;
    struct bed **flatten;
};

struct gtf_anno_type {
    enum exon_type type;
    // struct gtf *exon;
    int n, m;
    struct gene_type *a;
};

enum exon_type RE_type_map(uint8_t c);
const char *RE_tag_name(int i);
const char *exon_type_name(int i);
#endif
