#ifndef GTF_H
#define GTF_H

#include <stdlib.h>
#include "dict.h"

#define EXONIC    0
#define INTRONIC  1

enum feature_type {
    feature_unknow =-1,
    feature_gene,
    feature_transcript,
    feature_CDS,
    feature_start_codon,
    feature_stop_codon,
    feature_5UTR,
    feature_3UTR,
    feature_inter,
    feature_inter_CNS,
    feature_intron_CNS,
    feature_exon,
    feature_5UTR_alias,
    feature_3UTR_alias,
    feature_Selenocysteine,
};

struct gtf_lite {
    int seqname;
    int source;
    enum feature_type type;
    int start;
    int end;
    int strand; // 0 on forward, 1 on reverse
    int gene_id;
    int gene_name;
    int transcript_id;
    void *attr_dict;
    int n_son;
    int m_son;
    struct gtf_lite *son;
};

struct idx {
    int n, m;
    int *idx;
};
struct ctg_idx {
    int offset;
    int idx;
};
struct gtf_spec {
    struct dict *name; // contig names
    struct dict *gene_name;
    struct dict *gene_id;
    struct dict *transcript_id;
    struct dict *sources; //
    struct dict *attrs; // attributes
    struct dict *features;
    
    // index for contigs
    struct ctg_idx *ctg;
    uint64_t *idx;
    
    // index for genes
    // struct idx *gene_idx;
    
    int n_gtf, m_gtf;
    struct gtf_lite *gtf;    
};

const char *get_feature_name(enum feature_type type);
struct gtf_spec *gtf_spec_init();
void gtf_destory(struct gtf_spec *G);
struct gtf_spec *gtf_read(const char *fname);
struct gtf_lite *gtf_overlap_gene(struct gtf_spec *G, char *name, int start, int end, int *n, int cache);
/*
char *gtf_get_gene_name(struct gtf_spec *G, struct gtf_lite *gl);
char *gtf_get_gene_id(struct gtf_spec *G, struct gtf_lite *gl);
char *gtf_get_transcript_id(struct gtf_spec *G, struct gtf_lite *gl);
*/
#endif
