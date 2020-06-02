#ifndef GTF_H
#define GTF_H

#include <stdlib.h>
#include "dict.h"
#include "region_index.h"
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

struct gtf {
    int seqname;
    int source;
    enum feature_type type;
    int start;
    int end;
    int strand; // 0 on forward, 1 on reverse
    int gene_id;
    int gene_name;
    int transcript_id;
    struct dict *attr; // attributions
    struct dict *query; // used to fast access gtf
    int n_gtf, m_gtf;
    struct gtf *gtf;
};

struct _ctg_idx;

struct gtf_spec {
    struct dict *name; // contig names
    struct dict *gene_name;
    struct dict *gene_id;
    struct dict *transcript_id;
    struct dict *sources; //
    struct dict *attrs; // attributes
    struct dict *features;

    struct gtf_idx *idx;
    // index for contigs
    struct _ctg_idx *ctg;
    int n_gtf, m_gtf;
    struct gtf_lite *gtf;
};

struct gtf_ctg {
    struct dict *gene_idx;
    struct region_index *idx;    
    int n_gtf, m_gtf;
    struct gtf *gtf; 
};
struct gtf_spec2 {
    struct dict *name; // contig names
    struct dict *gene_name;
    struct dict *gene_id;
    struct dict *transcript_id;
    struct dict *sources; //
    struct dict *attrs; // attributes
    struct dict *features;
};

const char *get_feature_name(enum feature_type type);
struct gtf_spec *gtf_spec_init();
void gtf_destroy(struct gtf_spec *G);
struct gtf_spec *gtf_read(const char *fname, int filter);

struct region_itr *gtf_query(struct gtf_spec *G, char *name, int start, int end);

struct gtf_spec2 *gtf_read2(const char *fname, int filter);
struct region_itr *gtf_query2(struct gtf_spec2 *G, char *name, int start, int end);
void gtf_destroy2(struct gtf_spec2 *G);

#endif
