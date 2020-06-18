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

struct _ctg_idx {
    int offset;
    int idx;
};

const char *get_feature_name(enum feature_type type)
{
    assert(type>-1);
    return feature_type_names[type];
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
    int i;
    int l;
    l = sizeof(feature_type_names)/sizeof(feature_type_names[0]);
    for (i = 0; i < l; ++i) 
        dict_push(G->features, (char*)feature_type_names[i]);
        
    return G;
}
void gtf_reset(struct gtf *gtf)
{
    memset(gtf, 0, sizeof(struct gtf));
    gtf->seqname = gtf->source = gtf->start = gtf->end = gtf->gene_id = gtf->gene_name = gtf->transcript_id = -1;
}
static int cmpfunc (const void *_a, const void *_b)
{
    struct gtf_lite *a = (struct gtf_lite*)_a;
    struct gtf_lite *b = (struct gtf_lite*)_b;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    return a->end - b->end;
}
static int cmpfunc1 (const void *_a, const void *_b)
{
    struct gtf_lite *a = *(struct gtf_lite**)_a;
    struct gtf_lite *b = *(struct gtf_lite**)_b;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    return a->end - b->end;
}

struct gtf_idx {
    struct region_index *idx;
};

static void gtf_build_index(struct gtf_spec *G)
{
    // sort gene by coordinate
    qsort(G->gtf, G->n_gtf, sizeof(struct gtf_lite), cmpfunc);
    G->ctg = malloc(dict_size(G->name)*sizeof(struct _ctg_idx));
    memset(G->ctg, 0, sizeof(struct _ctg_idx)*dict_size(G->name));
    G->idx = malloc(dict_size(G->name)*sizeof(struct gtf_idx));    
    
    int i;
    for (i = 0; i < dict_size(G->name); ++i)
        G->idx[i].idx = region_index_build();
    
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *gl = &G->gtf[i];        
        G->ctg[gl->seqname].offset++;    
        if (G->ctg[gl->seqname].idx == 0) G->ctg[gl->seqname].idx = i+1;

        index_bin_push(G->idx[gl->seqname].idx, gl->start, gl->end, gl);
    }
    for (i = 0; i < dict_size(G->name); ++i) G->ctg[i].idx -= 1; // convert to 0 based

    // from v0.4, sort transcript and exon record in case GTF is not properly defined
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *g0 = &G->gtf[i]; // each gene
        qsort(g0->son, g0->n_son, sizeof(struct gtf_lite), cmpfunc);
        int j;
        for (j = 0; j < g0->n_son; ++j) { // each transcript
            struct gtf_lite *g1 = &g0->son[j];
            if (g1->n_son > 1) 
                qsort(g1->son, g1->n_son, sizeof(struct gtf_lite), cmpfunc);
        }
    }
}
static void gtf_lite_clean(struct gtf_lite *g)
{
    int i;
    for (i = 0; i < g->n_son; ++i) gtf_lite_clean(&g->son[i]);
    if (g->m_son) free(g->son);
    khint_t k;
    kh_attr_t *hash = (kh_attr_t*)g->attr_dict;
    for (k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (kh_exist(hash, k)) 
            if(kh_val(hash,k)) free(kh_val(hash,k));
    }
    kh_destroy(attr, hash);
}
void gtf_destroy(struct gtf_spec *G)
{
    int i;
    for (i = 0; i < G->n_gtf; ++i) gtf_lite_clean(&G->gtf[i]);

    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_idx *idx = &G->idx[i];
        region_index_destroy(idx->idx);
    }
    free(G->idx);
    free(G->ctg);
    free(G->gtf);
    dict_destroy(G->name);
    dict_destroy(G->gene_name);
    dict_destroy(G->gene_id);
    dict_destroy(G->transcript_id);
    dict_destroy(G->sources);
    dict_destroy(G->attrs);
    dict_destroy(G->features);
    free(G);
}

struct attr_pair {
    char *key;
    char *val;
};

static struct attr_pair *split_gff(kstring_t *str, int *_n)
{
    int i=0;
    int n =0, m = 0;
    struct attr_pair *pair = NULL;

    for (;;) {
        if (i >= str->l) break;
        if (n == m) {
            m += 4;
            pair = realloc(pair, sizeof(struct attr_pair)*m);
        }

        int j = str->l -1;
        while (isspace(str->s[j]) || str->s[j] == ';') j--;
        str->l = j+1;
        str->s[str->l] = '\0';
        
        kstring_t name = {0,0,0};
        kstring_t val = {0,0,0};
        while (i < str->l && !isspace(str->s[i]) && str->s[i] != ';') {
            kputc(str->s[i], &name);
            ++i;
        }

        while (i < str->l && isspace(str->s[i]) || str->s[i] == ';') ++i;
        
        if (str->s[i] == '"')  {
            ++i; // skip comma
            for (;i < str->l-1;) {
                if (str->s[i] == '"' && i +1 == str->l) {
                    i++; // skip ;
                    i++; // next record
                    break;
                }
                kputc(str->s[i], &val);
                i++;
            }
        }

        while (isspace(str->s[i])) ++i; // emit ends
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
    /*
    int *t = ksplit(&str, ';', n);
    struct attr_pair *pp = malloc(*n*sizeof(struct attr_pair));
    int i;
    kstring_t key = {0,0,0};
    for (i = 0; i < *n; ++i) {
        char *p0 = str.s+t[i];
        char *e = p0 + strlen(p0);
        while (isspace(*p0)) p0++;
        if (p0 == e || *p0 == '\0') break;
        char *p1 = p0;
        int j = 0;
        for (p1 = p0; !isspace(*p1) && p1 != e; p1++,j++);
        key.l = 0;
        kputsn(p0, j, &key);
        kputs("",&key);
        pp[i].key = strdup(key.s);
        if (*p1 == '\0') // flag
            pp[i].val = NULL;
        else {
            while (isspace(*p1)) p1++; // skip space
            p0 = p1;
            key.l = 0;
            if (*p0 != '"') { // no common
                for (j = 0, p1 = p0; *p1 != '\0'; ++p1,++j);
                kputsn(p0, j, &key);
                kputs("",&key);
                pp[i].val = strdup(key.s);
            }
            else {
                ++p0;
                for (j = 0, p1 = p0; *p1 != '"'; ++p1,++j);
                kputsn(p0, j, &key);
                kputs("",&key);
                pp[i].val = strdup(key.s);
            }
        }
    }
    *n = i; // in case empty endss
    free(str.s);
    free(key.s);
    free(t);
    return pp;
    */
}

static int gtf_push_new_gene(struct gtf_spec *G, struct gtf_lite *gl)
{
    if (G->n_gtf == G->m_gtf) {
        G->m_gtf = G->m_gtf == 0 ? 512 : G->m_gtf<<1;
        G->gtf = realloc(G->gtf, G->m_gtf*sizeof(struct gtf_lite));
    }
    // struct gtf_lite *g0 = &G->gtf[G->n_gtf];
    memcpy(&G->gtf[G->n_gtf], gl, sizeof(struct gtf_lite));
    G->n_gtf++;
    return 0;
}
static int gtf_push_to_record(struct gtf_lite *g0, struct gtf_lite *g1)
{
    if (g0->n_son == g0->m_son) {
        g0->m_son = g0->m_son == 0 ? 32 : g0->m_son<<1;
        g0->son = realloc(g0->son, g0->m_son*sizeof(struct gtf_lite));
    }
    memcpy(&g0->son[g0->n_son++], g1, sizeof(struct gtf_lite));
    return 0;
}
static int gtf_push_to_last_gene(struct gtf_spec *G, struct gtf_lite *gl)
{
    if (G->n_gtf == 0) error("No gene record found, bad format.");
    
    struct gtf_lite *g0 = &G->gtf[G->n_gtf-1];
    if (g0->type != feature_gene) error("Last record is not a gene, the GTF is not properly defined.");
    if (g0->gene_name != gl->gene_name) 
        error("g0->gene_name != gl->gene_name, %s vs %s", dict_name(G->gene_name, g0->gene_name), dict_name(G->gene_name, gl->gene_name));
    //assert(g0->gene_name == gl->gene_name);
    if (gl->type == feature_transcript) 
        gtf_push_to_record(g0, gl);
    else {
        if (g0->n_son == 0) error("No transcript record found, bad format");
        gtf_push_to_record(&g0->son[g0->n_son-1], gl);
    }
    return 0;
}

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
    
    if (filter && qry != feature_gene && qry != feature_exon && qry != feature_transcript && qry != feature_CDS) {
        free(s);
        return 0;
    }
    
    struct gtf_lite gtf;
    memset(&gtf, 0, sizeof(gtf));
    gtf.seqname = dict_push(G->name, str->s + s[0]);
    gtf.source = dict_push(G->sources, str->s + s[1]);
    gtf.type = qry;
    gtf.start = str2int(str->s+s[3]);
    gtf.end = str2int(str->s+s[4]);
    char *strand = str->s+s[6];
    gtf.strand = strand[0] == '+' ? 0 : 1;
    gtf.attr_dict = kh_init(attr);    
    char *attr = str->s+s[8];
    
    int i;
    int n0=0;
    struct attr_pair *pair = bend_pair(attr, &n0);
    for (i = 0; i < n0; ++i) {
        struct attr_pair *pp = &pair[i];
        if (strcmp(pp->key, "gene_id") == 0)
            gtf.gene_id = dict_push(G->gene_id, pp->val);       
        else if (strcmp(pp->key, "gene_name") == 0) 
            gtf.gene_name = dict_push(G->gene_name, pp->val);
        else if (strcmp(pp->key, "transcript_id") == 0) 
            gtf.transcript_id = dict_push(G->transcript_id, pp->val);
        else {
            int attr_id = dict_push(G->attrs, pp->key);
            khint_t k;
            int ret;
            k = kh_put(attr, (kh_attr_t*)gtf.attr_dict, attr_id, &ret);
            if (!ret) continue;
            kh_val((kh_attr_t*)gtf.attr_dict, k) = pp->val == NULL ? NULL : strdup(pp->val);
        }
        free(pp->key);
        if (pp->val) free(pp->val);
    }

    free(pair);
    free(s);

    // gtf_push(G, &gtf);    
    
    switch (qry) {
        case feature_gene:
        case feature_inter:
        case feature_inter_CNS:
            gtf_push_new_gene(G, &gtf);
            break;
        default:
            gtf_push_to_last_gene(G, &gtf);
            break;
    }
    
    return 0;
}

static kstring_t cache ={0,0,0};
static void reset_cache()
{
    cache.l = 0;
}
static void free_cache()
{
    free(cache.s);
}

static int gtf_push(struct gtf_spec2 *G, struct gtf_ctg *ctg, struct gtf *gtf, int feature)
{
    int gene_idx = dict_pushInt(ctg->gene_idx, gtf->gene_id);
    
    struct gtf *gene_gtf = dict_query_value(ctg->gene_idx, gene_idx);

    if (feature == feature_gene && gene_gtf != NULL) {
        warnings("Duplicated gene record? %s", dict_name(G->gene_name, gtf->gene_name));
        return 1;
    }
    
    if (gene_gtf == NULL) { // set new gene record
        if (ctg->n_gtf == ctg->m_gtf) {
            ctg->m_gtf = ctg->m_gtf == 0 ? 1000 : ctg->m_gtf*2;
            ctg->gtf = realloc(ctg->gtf, sizeof(struct gtf)*ctg->m_gtf);
        }
        
        gene_gtf = &ctg->gtf[ctg->n_gtf++];

        gtf_reset(gene_gtf);        

        gene_gtf->query = dict_init();
        dict_set_value(gene_gtf->query);
        dict_assign_value(ctg->gene_idx, gene_idx, gene_gtf);
        if (feature == feature_gene) {
            memcpy(gene_gtf, gtf, sizeof(struct gtf));
            return 0;
        }
    }

    if (gtf->transcript_id == -1) error("No transcript found.");
    
    // no gene record
    if (gene_gtf->query == NULL) {
        gene_gtf->query = dict_init();
        dict_set_value(gene_gtf->query);
    }

    if (gene_gtf->gene_id == -1) 
        gene_gtf->gene_id = gtf->gene_id;
    if (gene_gtf->gene_name == -1)
        gene_gtf->gene_name = gtf->gene_name;
    
    int trans_idx = dict_pushInt(gene_gtf->query, gtf->transcript_id);
    struct gtf *tx_gtf = dict_query_value(gene_gtf->query, trans_idx);

    if (feature == feature_transcript && tx_gtf != NULL) {
        warnings("Duplicated transcript record? %s", dict_name(G->transcript_id, gtf->transcript_id));
        return 1;
    }

    // setup transcript record
    if (tx_gtf == NULL) {
        if (gene_gtf->n_gtf == gene_gtf->m_gtf) {
            gene_gtf->m_gtf = gene_gtf->m_gtf == 0? 4 : gene_gtf->m_gtf*2;
            gene_gtf->gtf = realloc(gene_gtf->gtf, gene_gtf->m_gtf *sizeof(struct gtf));
        }
        tx_gtf = &gene_gtf->gtf[gene_gtf->n_gtf++];
        gtf_reset(tx_gtf);
        dict_assign_value(gene_gtf->query, trans_idx, tx_gtf);
        if (feature == feature_transcript) {
            memcpy(tx_gtf, gtf, sizeof(struct gtf));
            return 0; // trans record end here
        }
    }
    
    // exon, cds, UTRs etc.
    reset_cache();
    kputs(feature_type_names[feature], &cache); kputc(':', &cache);
    kputw(gtf->start, &cache); kputc('-', &cache);
    kputw(gtf->end, &cache);
    if (tx_gtf->query == NULL) {
        tx_gtf->query = dict_init();
        dict_set_value(tx_gtf->query);
    }
    int idx = dict_push(tx_gtf->query, cache.s);
    struct gtf *exon_gtf = dict_query_value(tx_gtf->query, idx);
    if (exon_gtf != NULL) {
        warnings("Duplicated record? %s", cache.s);
        return 1;
    }
    if (tx_gtf->n_gtf == tx_gtf->m_gtf) {
        tx_gtf->m_gtf = tx_gtf->m_gtf == 0 ? 4 : tx_gtf->m_gtf*2;
        tx_gtf->gtf = realloc(tx_gtf->gtf, tx_gtf->m_gtf*sizeof(struct gtf));
    }
    exon_gtf = &tx_gtf->gtf[tx_gtf->n_gtf++];
    memcpy(exon_gtf, gtf, sizeof(struct gtf));
    dict_assign_value(tx_gtf->query, idx, exon_gtf);

    return 0;
}

static int parse_str2(struct gtf_spec2 *G, kstring_t *str, int filter)
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
    
    if (filter &&
        qry != feature_gene &&
        qry != feature_exon &&
        qry != feature_transcript &&
        qry != feature_CDS &&
        qry != feature_5UTR &&
        qry != feature_3UTR) {
        free(s);
        return 0;
    }
    
    struct gtf gtf;
    gtf_reset(&gtf);
    gtf.seqname = dict_push(G->name, str->s + s[0]);
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
        else { // todo: update to dict structure
            int attr_id = dict_push(G->attrs, pp->key);
            if (gtf.attr == NULL) {
                gtf.attr = dict_init();
                dict_set_value(gtf.attr);
            }
            int idx = dict_pushInt(gtf.attr, attr_id);
            if (pp->val != NULL) {
                char *val = strdup(pp->val);
                dict_assign_value(gtf.attr, idx, val);
            }
        }
        free(pp->key);
        if (pp->val) free(pp->val);
    }
        
    free(pair);
    free(s);

    if (gtf.gene_id == -1 && gtf.gene_name == -1) {
        warnings("Record %s:%s:%d-%d has no gene_name and gene_id. Skip.", dict_name(G->name, gtf.seqname), feature_type_names[qry], gtf.start, gtf.end);
        return 1;
    }
    if (gtf.gene_id == -1) {
        warnings("Record %s:%s:%d-%d has no gene_id use gene_name instead.", dict_name(G->name, gtf.seqname), feature_type_names[qry], gtf.start, gtf.end);
        gtf.gene_id = dict_push(G->gene_id, dict_name(G->gene_name, gtf.gene_name));
    }

    gtf_push(G, ctg, &gtf, qry);

    return 0;
}
static void gtf_sort(struct gtf *gtf)
{
    int i;
    for (i = 0; i < gtf->n_gtf; ++i) 
        gtf_sort(&gtf->gtf[i]);
    if (gtf->n_gtf) {
        dict_destroy(gtf->query); // destroy query dict
        gtf->query = NULL;
        qsort((struct gtf*)gtf->gtf, gtf->n_gtf, sizeof(struct gtf), cmpfunc);
        int j;
        for (j = 0; j < gtf->n_gtf; ++j) {
            if (gtf->start > gtf->gtf[j].start) gtf->start = gtf->gtf[j].start;
            if (gtf->end < gtf->gtf[j].end) gtf->end = gtf->gtf[j].end;
        }
    }
}
static struct region_index *ctg_build_idx(struct gtf_ctg *ctg)
{
    struct region_index *idx = region_index_build();
    int i;
    for (i = 0; i < ctg->n_gtf; ++i) 
        index_bin_push(idx, ctg->gtf[i].start, ctg->gtf[i].end, &ctg->gtf[i]);
    return idx;
}
static int gtf_build_index2(struct gtf_spec2 *G)
{
    // update gene and transcript start and end record
    int i;
    int total_gene = 0;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name,i);
        assert(ctg);
        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            gtf_sort(&ctg->gtf[j]); // sort gene
            ctg->idx = ctg_build_idx(ctg);
        }
        total_gene+=ctg->n_gtf;
    }
    return total_gene;
}

// key names: gene_id, gene_name, transcript_id,
// if filter == 1, only keep genes, transcripts and exons
struct gtf_spec *gtf_read(const char *fname, int filter)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));

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
        parse_str(G, &str, filter); //warnings("Skip line %d, %s", line, str.s);

    }
    free(str.s);
    gzclose(fp);
    ks_destroy(ks);
    
    if (G->n_gtf == 0) {
        gtf_destroy(G);
        return NULL;
    }

    gtf_build_index(G);
    
    return G;
}

struct region_itr *gtf_query(struct gtf_spec *G, char *name, int start, int end)
{
    int id = dict_query(G->name, name);
    if (id == -1) return NULL;

    if (start < 0) start = 0;
    if (end < start) return NULL;

    int st = G->ctg[id].idx;
    if (end < G->gtf[st].start) return NULL; // out of range

    struct region_index *idx = G->idx[id].idx;

    struct region_itr *itr = region_query(idx, start, end);

    if (itr==NULL) return NULL;
    if (itr->n == 0) {
        free(itr);
        return NULL;
    }
    qsort((struct gtf_lite**)itr->rets, itr->n, sizeof(struct gtf_lite*), cmpfunc1);
    
    return itr;
}

struct gtf_spec2 *gtf_spec_init2()
{
    struct gtf_spec2 *G = malloc(sizeof(*G));
    memset(G, 0, sizeof(*G));
    G->name            = dict_init();
    G->gene_name       = dict_init();
    G->gene_id         = dict_init();
    G->transcript_id   = dict_init();
    G->sources         = dict_init();
    G->attrs           = dict_init();
    G->features        = dict_init();

    dict_set_value(G->name);
    int i;
    int l;
    l = sizeof(feature_type_names)/sizeof(feature_type_names[0]);
    for (i = 0; i < l; ++i) 
        dict_push(G->features, (char*)feature_type_names[i]);
    
    return G;
}
struct gtf_spec2 *gtf_read2(const char *fname, int f)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int line = 0;
    struct gtf_spec2 *G = gtf_spec_init2();
    
    while (ks_getuntil(ks, 2, &str, &ret)>=0) {
        line++;
        if (str.l == 0) {
            warnings("Line %d is empty. Skip.", line);
            continue;
        }
        if (str.s[0] == '#') continue;
        parse_str2(G, &str, f); //warnings("Skip line %d, %s", line, str.s);
    }
    free(str.s);
    gzclose(fp);
    ks_destroy(ks);
    
    if (dict_size(G->name) == 0) {
        gtf_destroy2(G);
        return NULL;
    }

    int n_gene = gtf_build_index2(G);
    LOG_print("Load %d genes.", n_gene);
    free_cache();
    return G;

}

struct region_itr *gtf_query2(struct gtf_spec2 const *G, char *name, int start, int end)
{
    int id = dict_query(G->name, name);
    if (id == -1) return NULL;

    if (start < 0) start = 0;
    if (end < start) return NULL;

    struct gtf_ctg *ctg = dict_query_value(G->name, id);
    if (ctg->n_gtf == 0) return NULL; // empty, should not happen?
    if (end < ctg->gtf[0].start) return NULL; // out of range

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
void gtf_clear(struct gtf *gtf)
{
    int i;
    for (i = 0; i < gtf->n_gtf; ++i)
        gtf_clear(&gtf->gtf[i]);
    if (gtf->n_gtf) free(gtf->gtf);

    if (gtf->attr != NULL) {
        int i;
        for (i = 0; i < dict_size(gtf->attr); ++i) {
            char *val = dict_query_value(gtf->attr, i);
            if (val) free(val);
        }
        dict_destroy(gtf->attr);
    }

    if (gtf->query) // usually already freed during indexing
        dict_destroy(gtf->query);
}
void gtf_destroy2(struct gtf_spec2 *G)
{
    int i;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        if (ctg->gene_idx) dict_destroy(ctg->gene_idx);
        region_index_destroy(ctg->idx);
        
        int j;
        for (j = 0; j < ctg->n_gtf; ++j)
            gtf_clear(&ctg->gtf[j]);
        
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

#ifdef GTF_MAIN
int main(int argc, char **argv)
{
    if (argc != 2) error("gtfformat in.gtf");
    struct gtf_spec *G = gtf_read(argv[1], 1);
    gtf_format_print_test(G);
    gtf_destory(G);
    return 0;
}

#endif
