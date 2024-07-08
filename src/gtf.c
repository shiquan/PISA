#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/ksort.h"
#include "dict.h"
#include "gtf.h"
#include "region_index.h"
#include "number.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8193)

KHASH_MAP_INIT_INT(attr, char*)

struct _ctg_idx {
    int offset;
    int idx;
};

static const char *feature_type_names[] = {
    // The following feature types are required: "gene", "transcript"
    "gene",
    "transcript",
    // The features "CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional.
    "CDS",
    "start_codon",
    "stop_codon",
    "5UTR",
    "3UTR",
    "inter",
    "inter_CNS",
    "intron_CNS",
    "exon",
    "five_prime_utr",
    "three_prime_utr",
    "Selenocysteine"
    // All other features will be ignored. The types must have the correct capitalization shown here.
};

const char *get_feature_name(enum feature_type type)
{
    assert(type>-1);
    return feature_type_names[type];
}

void gtf_reset(struct gtf *gtf)
{
    memset(gtf, 0, sizeof(struct gtf));
    gtf->seqname = gtf->source = gtf->start = gtf->end = gtf->gene_id = gtf->gene_name = gtf->transcript_id = -1;
}
struct gtf *gtf_create()
{
    struct gtf *g = malloc(sizeof(*g));
    gtf_reset(g);
    return g;
}
void gtf_attr_clear(struct attr *attr)
{
    if (attr) {
        if (attr->next) gtf_attr_clear(attr->next);
        if (attr->val) free(attr->val);
        free(attr);
    }
}
void gtf_clear(struct gtf *gtf)
{
    for (int i = 0; i < gtf->n_gtf; ++i) {
        gtf_clear(gtf->gtf[i]);
        free(gtf->gtf[i]);
    }
    if (gtf->n_gtf) free(gtf->gtf);
    if (gtf->attr != NULL) gtf_attr_clear(gtf->attr);
}

void gtf_copy(struct gtf *dest, struct gtf *src)
{
    gtf_reset(dest);
    dest->seqname = src->seqname;
    dest->source = src->source;
    dest->type = src->type;
    dest->start = src->start;
    dest->end = src->end;
    dest->strand = src->strand;
    dest->gene_id = src->gene_id;
    dest->gene_name = src->gene_name;
    dest->transcript_id = src->transcript_id;
    // copy attributes
    dest->attr = src->attr;

    src->attr = NULL; // do we need transfer the point??
}
static int cmpfunc1 (const void *_a, const void *_b)
{
    const struct gtf *a = *(const struct gtf**)_a;
    const struct gtf *b = *(const struct gtf**)_b;
    if (a->seqname != b->seqname) return (a->seqname > b->seqname) - (a->seqname < b->seqname);
    if (a->start != b->start) return (a->start > b->start) - (a->start < b->start);
    return (a->end > b->end) - (a->end < b->end) ;
}

struct gtf_idx {
    struct region_index *idx;
};

struct attr_pair {
    char *key;
    char *val;
};

static struct attr_pair *split_gff(kstring_t *str, int *_n)
{
    int i=0;
    int n =0, m = 0;
    struct attr_pair *pair = NULL;

    int j = str->l -1;
    while (isspace(str->s[j]) || str->s[j] == ';') j--;
    str->l = j+1;
    str->s[str->l] = '\0';

    for (;;) {
        if (i >= str->l) break;
        if (n == m) {
            m += 4;
            pair = realloc(pair, sizeof(struct attr_pair)*m);
        }
        
        kstring_t name = {0,0,0};
        kstring_t val = {0,0,0};
        while (i < str->l && !isspace(str->s[i]) && str->s[i] != ';') {
            kputc(str->s[i], &name);
            ++i;
        }

        while (isspace(str->s[i]) || str->s[i] == ';') ++i; // emit middle spaces
        
        if (str->s[i] == '"')  {
            ++i; // skip comma
            for (;i < str->l;) {
                if (str->s[i] == '"' || i+1 == str->l) {
                    i++; // skip ;
                    i++; // next record
                    break;
                }
                kputc(str->s[i], &val);
                i++;
            }
        }

        while (i < str->l && (isspace(str->s[i]) || str->s[i] == ';')) ++i; // emit ends
        if (name.l == 0) {
            warnings("Empty key. %s", str->s);
            continue;
        }

        pair[n].key = name.s;
        pair[n].val = val.s;
        n++;
    }
    *_n = n;
    return pair;
}
static struct attr_pair *bend_pair(char *s, int *n)
{    
    if (s == NULL) return NULL;
    kstring_t str = {0,0,0};    
    kputs(s, &str);

    struct attr_pair *p = split_gff(&str, n);
    free(str.s);
    return p;
}
static int gtf_push(struct gtf_spec *G, struct gtf_ctg *ctg, struct gtf *gtf, int feature)
{
    char *gene_id = dict_name(G->gene_id, gtf->gene_id);
    int gene_idx = dict_push(ctg->gene_idx, gene_id);
        
    struct gtf *gene_gtf = dict_query_value(ctg->gene_idx, gene_idx);
    //struct gtf *gene_gtf = dict_query_value(G->gene_id, gtf->gene_id);

    if (feature == feature_gene && gene_gtf != NULL) {
        warnings("Duplicated gene record? %s", dict_name(G->gene_name, gtf->gene_name));
        if (gene_gtf->start < 0 || gene_gtf->start > gtf->start) gene_gtf->start = gtf->start;
        if (gene_gtf->end < 0 || gene_gtf->end < gtf->end) gene_gtf->end = gtf->end;

        gene_gtf->strand = gtf->strand;
        gene_gtf->gene_id = gtf->gene_id;
        gene_gtf->gene_name = gtf->gene_name;
        gene_gtf->seqname = gtf->seqname;
        return 1;
    }
    
    if (gene_gtf == NULL) { // set new gene record
        if (ctg->n_gtf == ctg->m_gtf) {
            ctg->m_gtf = ctg->m_gtf == 0 ? 4: ctg->m_gtf*2;
            ctg->gtf = realloc(ctg->gtf, sizeof(struct gtf*)*ctg->m_gtf);
        }
        ctg->gtf[ctg->n_gtf] = gtf_create();
        gene_gtf = ctg->gtf[ctg->n_gtf++];

        // assign to ctg, because will reorder these records during indexing
        dict_assign_value(ctg->gene_idx, gene_idx, gene_gtf);
        //dict_assign_value(G->gene_id, gtf->gene_id, gene_gtf);
        //dict_assign_value(G->gene_name, gtf->gene_name, gene_gtf);
        
        gtf_reset(gene_gtf); 
        if (feature == feature_gene) {
            //memcpy(gene_gtf, gtf, sizeof(struct gtf));
            gtf_copy(gene_gtf, gtf);
            return 0;
        }
        gene_gtf->type = feature_gene;
        gene_gtf->strand = gtf->strand;
        gene_gtf->seqname = gtf->seqname;
    }

    if (gene_gtf->start < 0 || gene_gtf->start > gtf->start) gene_gtf->start = gtf->start;
    if (gene_gtf->end < 0 || gene_gtf->end < gtf->end) gene_gtf->end = gtf->end;

    if (gtf->transcript_id == -1) error("No transcript found. %s:%s:%d:%d", feature_type_names[feature], dict_name(G->name, gtf->seqname), gtf->start, gtf->end);

    // no gene record
    /*
    if (gene_gtf->query == NULL) {
        gene_gtf->query = dict_init();
        dict_set_value(gene_gtf->query);
    }
    */
    // inhert gene name
    if (gene_gtf->gene_id == -1) gene_gtf->gene_id = gtf->gene_id;
    if (gene_gtf->gene_name == -1) gene_gtf->gene_name = gtf->gene_name;

    //char *transcript_id = dict_query(G->transcript_id, gtf->transcript_id);
    //int trans_idx = dict_push(gene_gtf->query, gtf->transcript_id);
    //struct gtf *tx_gtf = dict_query_value(gene_gtf->query, trans_idx);
    struct gtf *tx_gtf = NULL; // dict_query_value(G->transcript_id, gtf->transcript_id);
    int i;
    for (i = 0; i < gene_gtf->n_gtf; ++i) {
        struct gtf *tx_gtf0 = gene_gtf->gtf[i];
        if (gtf->transcript_id == tx_gtf0->transcript_id) {
            tx_gtf = tx_gtf0;
            break;
        }
    }
    if (feature == feature_transcript && tx_gtf != NULL) {
        warnings("Duplicated transcript record? %s", dict_name(G->transcript_id, gtf->transcript_id));
        if (tx_gtf->start < 0 || tx_gtf->start > gtf->start) tx_gtf->start = gtf->start;
        if (tx_gtf->end < 0 || tx_gtf->end < gtf->end) tx_gtf->end = gtf->end;
        tx_gtf->strand = gtf->strand;
        tx_gtf->gene_id = gtf->gene_id;
        tx_gtf->gene_name = gtf->gene_name;
        tx_gtf->seqname = gtf->seqname;
        return 1;
    }

    // setup transcript record
    if (tx_gtf == NULL) {
        if (gene_gtf->n_gtf == gene_gtf->m_gtf) {
            gene_gtf->m_gtf = gene_gtf->m_gtf == 0? 4 : gene_gtf->m_gtf*2;
            gene_gtf->gtf = realloc(gene_gtf->gtf, gene_gtf->m_gtf *sizeof(struct gtf*));
        }
        gene_gtf->gtf[gene_gtf->n_gtf] = gtf_create();
        tx_gtf = gene_gtf->gtf[gene_gtf->n_gtf++];
        //dict_assign_value(gene_gtf->query, trans_idx,tx_gtf);
        // dict_assign_value(G->transcript_id, gtf->transcript_id, tx_gtf);
        
        gtf_reset(tx_gtf);

        if (feature == feature_transcript) {
            //memcpy(tx_gtf, gtf, sizeof(struct gtf));
            gtf_copy(tx_gtf, gtf);
            
            return 0; // trans record end here
        }

        // no transcript record in the gtf, produce a record automatically
        // init gtf type
        tx_gtf->type = feature_transcript;
        tx_gtf->strand = gtf->strand;
        tx_gtf->seqname = gtf->seqname;        
    }

    // inhert gene and transcript name 
    if (tx_gtf->gene_id == -1) tx_gtf->gene_id = gtf->gene_id;
    if (tx_gtf->gene_name == -1) tx_gtf->gene_name = gtf->gene_name;
    if (tx_gtf->transcript_id == -1) tx_gtf->transcript_id = gtf->transcript_id;

    // will update before indexing
    if (tx_gtf->start < 0 || tx_gtf->start > gtf->start) tx_gtf->start = gtf->start;
    if (tx_gtf->end < 0 || tx_gtf->end < gtf->end) tx_gtf->end = gtf->end;
    
    if (tx_gtf->gene_id != gtf->gene_id) {
        error("Inconsitance gene id between transcript %s and gene %s.",
              dict_name(G->transcript_id, tx_gtf->transcript_id),
              dict_name(G->gene_id, gtf->gene_id)
            );
    }
    if (tx_gtf->gene_name != gtf->gene_name) {
        error("Inconsitance gene name between transcript %s and gene %s.",
              dict_name(G->transcript_id, tx_gtf->transcript_id),
              dict_name(G->gene_name, gtf->gene_name)
            );        
    }
    if (tx_gtf->transcript_id != gtf->transcript_id) {
        error("Inconsitance transcript id between transcript %s and record %s.",
              dict_name(G->transcript_id, tx_gtf->transcript_id),
              dict_name(G->transcript_id, gtf->transcript_id)
            );        
    }
    // assert(tx_gtf->gene_id == gtf->gene_id);
    // assert(tx_gtf->gene_name == gtf->gene_name);
    // assert(tx_gtf->transcript_id == gtf->transcript_id);

    // exon, cds, UTRs etc.
    /*
      reset_cache();
    kputs(feature_type_names[feature], &cache); kputc(':', &cache);
    kputw(gtf->start, &cache); kputc('-', &cache);
    kputw(gtf->end, &cache);
    if (tx_gtf->query == NULL) {
        tx_gtf->query = dict_init();
        dict_set_value(tx_gtf->query);
    }
    //int idx = dict_push(tx_gtf->query, cache.s);
    struct gtf *exon_gtf = dict_query_value(tx_gtf->query, idx);
    if (exon_gtf != NULL) {
        warnings("Duplicated record? %s", cache.s);
        return 1;
    }
    */
    if (tx_gtf->n_gtf == tx_gtf->m_gtf) {
        tx_gtf->m_gtf = tx_gtf->m_gtf == 0 ? 4 : tx_gtf->m_gtf*2;
        tx_gtf->gtf = realloc(tx_gtf->gtf, tx_gtf->m_gtf*sizeof(struct gtf*));
    }
    tx_gtf->gtf[tx_gtf->n_gtf] = gtf_create();
    struct gtf *exon_gtf = tx_gtf->gtf[tx_gtf->n_gtf++];
    //memcpy(exon_gtf, gtf, sizeof(struct gtf));
    
    gtf_copy(exon_gtf, gtf);
    //dict_assign_value(tx_gtf->query, idx, exon_gtf);


    // update coding signature
    if (exon_gtf->type == feature_CDS) {
        tx_gtf->coding = 1;
        gene_gtf->coding = 1;
    }
    
    return 0;
}

#define FILTER_ATTRS  2
#define FILTER_TRANS  1

static int parse_str(struct gtf_spec *G, kstring_t *str, int filter)
{
    int n;
    int *s = ksplit(str, '\t', &n);
    if (n != 9) error("Unknown format. %s", str->s);
    
    char *feature = str->s + s[2];

    int qry = dict_query(G->features, feature);
    if (qry == -1) {
        free(s);
        return 1;
    }
    
    if ((filter &0x3) & FILTER_TRANS) {
        if (qry != feature_gene && qry != feature_exon &&
            qry != feature_transcript) {
            //qry != feature_CDS &&
            //qry != feature_5UTR &&
            //qry != feature_3UTR) {
            free(s);
            return 0;
        }
    }

    struct gtf gtf;
    gtf_reset(&gtf);
    gtf.seqname = dict_push(G->name, str->s + s[0]);
    if (s[1] != 1 || str->s[s[1]] != '.')
        gtf.source = dict_push(G->sources, str->s + s[1]);
    gtf.type = qry;
    gtf.start = str2int(str->s+s[3]);
    gtf.end = str2int(str->s+s[4]);
    char *strand = str->s+s[6];
    gtf.strand = strand[0] == '-' ? 1 : 0;
    char *attr = str->s+s[8];

    struct gtf_ctg *ctg = dict_query_value(G->name, gtf.seqname);
    if (ctg == NULL) { // init contig value
        ctg = malloc(sizeof(struct gtf_ctg));
        memset(ctg, 0, sizeof(struct gtf_ctg));
        ctg->gene_idx = dict_init();
        dict_set_value(ctg->gene_idx);
        dict_assign_value(G->name, gtf.seqname, ctg);
    }

    int i;
    int n0=0;
    struct attr_pair *pair = bend_pair(attr, &n0);
    for (i = 0; i < n0; ++i) {
        struct attr_pair *pp = &pair[i];
        if (strcmp(pp->key, "gene_id") == 0)
            gtf.gene_id = dict_push(G->gene_id, pp->val);       
        else if (strcmp(pp->key, "gene_name") == 0)
            gtf.gene_name = dict_push(G->gene_name, pp->val);
        else if (strcmp(pp->key, "gene") == 0) // some gtf use gene instead of gene_name
            gtf.gene_name = dict_push(G->gene_name, pp->val);
        else if (strcmp(pp->key, "transcript_id") == 0) 
            gtf.transcript_id = dict_push(G->transcript_id, pp->val);
        else {
            if ((filter & 0x3) & FILTER_ATTRS) { // todo: update to dict structure
                //int attr_id = dict_push(G->attrs, pp->key);
                /* dict_push(G->attrs, pp->key); */
                /* if (gtf.attr == NULL) { */
                /*     gtf.attr = dict_init(); */
                /*     dict_set_value(gtf.attr); */
                /* } */
                /* int idx = dict_push(gtf.attr, pp->key); */
                /* if (pp->val != NULL) { */
                /*     char *val = strdup(pp->val); */
                /*     dict_assign_value(gtf.attr, idx, val); */
                /* } */
                                
                struct attr *attr = malloc(sizeof(struct attr));
                attr->id = dict_push(G->attrs, pp->key);
                attr->val = NULL;
                attr->next = NULL;
                if (gtf.attr == NULL) gtf.attr = attr;
                else {
                    struct attr *tmp;
                    for (tmp = gtf.attr; tmp->next; tmp = tmp->next);
                    tmp->next = attr;
                }

                if (pp->val != NULL) attr->val = strdup(pp->val);
            }
        }
        free(pp->key);
        if (pp->val) free(pp->val);
    }
        
    free(pair);
    free(s);

    if (gtf.gene_id == -1 && gtf.gene_name == -1) {
        warnings("Record %s:%s:%d-%d has no gene_name and gene_id. Skip.", dict_name(G->name, gtf.seqname), feature_type_names[qry], gtf.start, gtf.end);
        gtf_clear(&gtf);
        return 1;
    }
    if (gtf.gene_id == -1) {
        // warnings("Record %s:%s:%d-%d has no gene_id, use gene_name instead.", dict_name(G->name, gtf.seqname), feature_type_names[qry], gtf.start, gtf.end);
        gtf.gene_id = dict_push(G->gene_id, dict_name(G->gene_name, gtf.gene_name));
    }

    if (gtf.gene_name == -1) {
        // warnings("Record %s:%s:%d-%d has no gene_name, use gene_id instead.", dict_name(G->name, gtf.seqname), feature_type_names[qry], gtf.start, gtf.end);
        gtf.gene_name = dict_push(G->gene_name, dict_name(G->gene_id, gtf.gene_id));
    }
    if (gtf_push(G, ctg, &gtf, qry)) {
        warnings("Failed to push record, %s:%s:%d-%d", dict_name(G->name, gtf.seqname), feature_type_names[qry], gtf.start, gtf.end);
    }
    gtf_clear(&gtf);
    return 0;
}
static void gtf_sort(struct gtf *gtf)
{
    int i;
    for (i = 0; i < gtf->n_gtf; ++i)
        gtf_sort(gtf->gtf[i]);
        
    if (gtf->n_gtf) {
        qsort((const struct gtf**)gtf->gtf, gtf->n_gtf, sizeof(struct gtf*), cmpfunc1);
        assert(gtf->start < gtf->end);
    }
    /*
    if (gtf->query) {
        dict_destroy(gtf->query); // destroy query dict
        gtf->query = NULL;
    }
    */

}
static struct region_index *ctg_build_idx(struct gtf_ctg *ctg)
{
    struct region_index *idx = region_index_create();
    int i;
    int last_start = -1;
    int last_id = -1;
        
    for (i = 0; i < ctg->n_gtf; ++i) {
        if (last_id != ctg->gtf[i]->seqname) {
            last_id = ctg->gtf[i]->seqname;
            last_start = -1;
        } else {
            if (last_start > ctg->gtf[i]->start) {
                error("gtf is not sorted.");   
            }
        }
        last_start = ctg->gtf[i]->start;
        index_bin_push(idx, ctg->gtf[i]->start, ctg->gtf[i]->end, ctg->gtf[i]);
    }
    return idx;
}
static int gtf_build_index(struct gtf_spec *G)
{
    // update gene and transcript start and end record
    int i;
    int total_gene = 0;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name,i);
        assert(ctg);
        qsort((const struct gtf**)ctg->gtf, ctg->n_gtf, sizeof(struct gtf*), cmpfunc1);
        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            gtf_sort(ctg->gtf[j]); // sort gene
            struct gtf *g = ctg->gtf[j];
            struct gtf *g0 = dict_query_value(G->gene_name, g->gene_name);
            if (g0 == NULL) {
                dict_assign_value(G->gene_name, g->gene_name, g);
            } else {
                for (;;) {
                    if (g0->ext == NULL) {
                        g0->ext = g;
                        break;
                    }
                    g0 = g0->ext;
                }
            }
            
            for (int k = 0; k < g->n_gtf; ++k) {
                struct gtf *tx = g->gtf[k];
                struct gtf *tx0 = dict_query_value(G->transcript_id, tx->transcript_id);
                if (tx0 == NULL) {
                    dict_assign_value(G->transcript_id, tx->transcript_id, tx);
                } else {
                    for (;;) {
                        if (tx0->ext == NULL) {
                            tx0->ext = tx;
                            break;
                        }
                        tx0 = tx0->ext;
                    }
                }
            }
        }
        ctg->idx = ctg_build_idx(ctg);
        total_gene+=ctg->n_gtf;
    }
    return total_gene;
}

struct gtf_spec *gtf_spec_init()
{
    struct gtf_spec *G = malloc(sizeof(*G));
    memset(G, 0, sizeof(*G));
    G->name            = dict_init();
    G->gene_name       = dict_init();
    G->gene_id         = dict_init();
    G->transcript_id   = dict_init();
    G->sources         = dict_init();
    G->attrs           = dict_init();
    G->features        = dict_init();

    dict_set_value(G->name);
    dict_set_value(G->gene_name);
    //dict_set_value(G->gene_id);
    dict_set_value(G->transcript_id);
    int i;
    int l;
    l = sizeof(feature_type_names)/sizeof(feature_type_names[0]);
    for (i = 0; i < l; ++i) 
        dict_push(G->features, (char*)feature_type_names[i]);

    // todo: more features..
    
    return G;
}

struct gtf_spec *gtf_read(const char *fname, int f)
{
    LOG_print("GTF loading..");
    double t_real;
    t_real = realtime();

    gzFile fp;
    fp = gzopen(fname, "r");
    if(fp == NULL) error("%s : %s.", fname, strerror(errno));

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int line = 0;
    struct gtf_spec *G = gtf_spec_init();
    
    while (ks_getuntil(ks, 2, &str, &ret)>=0) {
        line++;
        if (str.l == 0) {
            warnings("Line %d is empty. Skip.", line);
            continue;
        }
        if (str.s[0] == '#') continue;
        parse_str(G, &str, f); //warnings("Skip line %d, %s", line, str.s);
    }
    free(str.s);
    gzclose(fp);
    ks_destroy(ks);
    
    if (dict_size(G->name) == 0) {
        gtf_destroy(G);
        return NULL;
    }

    int n_gene = gtf_build_index(G);
    LOG_print("Load %d genes.", n_gene);
    //free_cache();
    LOG_print("Load time : %.3f sec", realtime() - t_real);
    return G;

}

struct gtf_spec *gtf_read_lite(const char *fname)
{
    return gtf_read(fname, 1);
}
struct region_itr *gtf_query(struct gtf_spec const *G, const char *name, int start, int end)
{
    int id = dict_query(G->name, name);
    if (id == -1) return NULL;

    if (start < 0) start = 0;
    // if (end < start) return NULL;
    if (end < -1) return NULL;
    if (end == -1) end = INT_MAX;
    
    struct gtf_ctg *ctg = dict_query_value(G->name, id);
    if (ctg->n_gtf == 0) return NULL; // empty, should not happen?
    if (end < ctg->gtf[0]->start) return NULL; // out of range
    
    struct region_index *idx = ctg->idx;

    struct region_itr *itr = region_query(idx, start, end);

    if (itr==NULL) return NULL;
    if (itr->n == 0) {
        free(itr);
        return NULL;
    }
    qsort((struct gtf**)itr->rets, itr->n, sizeof(struct gtf*), cmpfunc1);
    
    return itr;
}
void gtf_destroy(struct gtf_spec *G)
{
    int i;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        if (ctg->gene_idx) dict_destroy(ctg->gene_idx);
        region_index_destroy(ctg->idx);
        
        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            gtf_clear(ctg->gtf[j]);
            free(ctg->gtf[j]);
        }
        
        free(ctg->gtf);
        free(ctg);
    }
    dict_destroy(G->name);
    dict_destroy(G->gene_name);
    dict_destroy(G->gene_id);
    dict_destroy(G->transcript_id);
    dict_destroy(G->sources);
    dict_destroy(G->attrs);
    dict_destroy(G->features);
    free(G);
}

char *GTF_seqname(struct gtf_spec *G, int id)
{
    return dict_name(G->name, id);
}

char *GTF_genename(struct gtf_spec *G, int id)
{
    return dict_name(G->gene_name, id);
}
char *GTF_transid(struct gtf_spec *G, int id)
{
    return dict_name(G->transcript_id, id);
}

void write_gtf_fp(struct gtf_spec *G, struct gtf *gtf, FILE *fp, struct dict *keys)
{
    fputs(dict_name(G->name,gtf->seqname), fp);
    if (gtf->source == -1)fputs("\t.\t", fp);
    else fprintf(fp, "\t%s\t", dict_name(G->sources, gtf->source));
    fputs(feature_type_names[gtf->type], fp);
    fprintf(fp,"\t%d\t%d\t.\t%c\t.\t", gtf->start, gtf->end, "+-"[gtf->strand]);
    fprintf(fp, "gene_name \"%s\"; gene_id \"%s\";", dict_name(G->gene_name, gtf->gene_name),
            dict_name(G->gene_id, gtf->gene_id));

    if (gtf->transcript_id >= 0) {
        fprintf(fp, " transcript_id \"%s\";", dict_name(G->transcript_id, gtf->transcript_id));
    }
    
    if (gtf->attr) {
        struct attr *tmp = gtf->attr;
        for (tmp = gtf->attr; tmp; tmp = tmp->next) {
            char *name = dict_name(G->attrs, tmp->id);
            fprintf(fp, " %s", name);
            if (tmp->val) fprintf(fp, " \"%s\";", tmp->val);
            else fputc(';', fp);
        }
    }

    fputc('\n', fp);
    
    int i;
    for (i = 0; i < gtf->n_gtf; ++i) write_gtf_fp(G, gtf->gtf[i], fp, keys);
}

void gtf_dump(struct gtf_spec *G, const char *fname, struct dict *keys)
{
    FILE *fp = fname == NULL ? stdout : fopen(fname, "w");
    if (fp == NULL) error("%s : %s.", fname, strerror(errno));
    int i, j;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gtf = ctg->gtf[j];
            write_gtf_fp(G, gtf, fp, keys);
        }
    }
    if (fp != stdout) fclose(fp);
}

struct gtf *gtf_query_gene(struct gtf_spec *G, const char *name)
{
    return (struct gtf*)dict_query_value2(G->gene_name, name);
}

struct gtf *gtf_query_tx(struct gtf_spec *G, const char *name)
{
    return (struct gtf*)dict_query_value2(G->transcript_id, name);
}

#ifdef GTF_MAIN
int main(int argc, char **argv)
{
    if (argc != 2) error("gtfformat in.gtf");
    struct gtf_spec *G = gtf_read_lite(argv[1]);
    //gtf_format_print_test(G);
    gtf_destroy(G);
    return 0;
}

#endif
