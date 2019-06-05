#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/ksort.h"
#include "dict.h"
#include "gtf.h"
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
static int cmpfunc (const void *_a, const void *_b)
{
    struct gtf_lite *a = (struct gtf_lite*)_a;
    struct gtf_lite *b = (struct gtf_lite*)_b;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    return a->end - b->end;
}

static void gtf_build_index(struct gtf_spec *G)
{
    qsort(G->gtf, G->n_gtf, sizeof(struct gtf_lite), cmpfunc);
    G->ctg = malloc(G->name->n*sizeof(struct ctg_idx));
    memset(G->ctg, 0, sizeof(struct ctg_idx)*G->name->n);
    G->idx = malloc(G->n_gtf*sizeof(uint64_t));
    int i, j;

    for (i = 0, j = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *gl = &G->gtf[i];
        G->idx[i] = (uint64_t)gl->start<<32|gl->end;
        G->ctg[gl->seqname].offset++;        
        if (G->ctg[gl->seqname].idx == 0) G->ctg[gl->seqname].idx = i+1;
    }
    for (i = 0; i < G->n_gtf; ++i) G->ctg[i].idx -= 1; // convert to 0 based
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
void gtf_destory(struct gtf_spec *G)
{
    int i;
    for (i = 0; i < G->n_gtf; ++i) gtf_lite_clean(&G->gtf[i]);
    free(G->ctg);
    free(G->idx);
    // free(G->gene_idx);
    free(G->gtf);
    dict_destory(G->name);
    dict_destory(G->gene_name);
    dict_destory(G->gene_id);
    dict_destory(G->transcript_id);
    dict_destory(G->sources);
    dict_destory(G->attrs);
    dict_destory(G->features);
    free(G);
}

struct attr_pair {
    char *key;
    char *val;
};
static struct attr_pair *bend_pair(char *s, int *n)
{
    if (s == NULL) return NULL;
    kstring_t str = {0,0,0};    
    kputs(s, &str);
    int *t = ksplit(&str, ';', n);
    struct attr_pair *pp = malloc(*n*sizeof(struct attr_pair));
    int i;
    kstring_t key = {0,0,0};
    for (i = 0; i < *n; ++i) {
        char *p0 = str.s+t[i];
        while (isspace(*p0)) p0++;
        char *p1 = p0;
        int j = 0;
        for (p1 = p0; !isspace(*p1); p1++,j++);
        key.l = 0;
        kputsn(p0, j, &key);
        kputs("",&key);
        pp[i].key = strdup(key.s);
        if (*p1 == '\0') // flag
            pp[i].val = NULL;
        else {
            p0 = ++p1;
            key.l = 0;
            if (*p0 != '"') error("Bad format.");
            ++p0;
            for (j = 0, p1 = p0; *p1 != '"'; ++p1,++j);
            kputsn(p0, j, &key);
            kputs("",&key);
            pp[i].val = strdup(key.s);
        }        
    }
    free(str.s);
    free(key.s);
    free(t);
    return pp;
}

static int gtf_push_new_gene(struct gtf_spec *G, struct gtf_lite *gl)
{
    if (G->n_gtf == G->m_gtf) {
        G->m_gtf = G->m_gtf == 0 ? 512 : G->m_gtf<<1;
        G->gtf = realloc(G->gtf, G->m_gtf*sizeof(struct gtf_lite));
    }
    memcpy(&G->gtf[G->n_gtf], gl, sizeof(struct gtf_lite));
    G->n_gtf++;
    return 0;
}
static int gtf_push_to_record(struct gtf_lite *g0, struct gtf_lite *g1)
{
    if (g0->n_son == g0->m_son) {
        g0->m_son = g0->m_son == 0 ? 2 : g0->m_son<<1;
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
    if (gl->type == feature_transcript) 
        gtf_push_to_record(g0, gl);
    else {
        if (g0->n_son == 0) error("No transcript record found, bad format");
        gtf_push_to_record(&g0->son[g0->n_son-1], gl);
    }
    return 0;
}

/*
1	ensembl_havana	gene	3205901	3671498	.	-	.	gene_id "ENSMUSG00000051951"; gene_version "5"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"
1	havana	transcript	3205901	3216344	.	-	.	gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "Xkr4-203"; transcript_source "havana"; transcript_biotype "processed_transcript"; transcript_support_level "1"
 */
static int parse_str(struct gtf_spec *G, kstring_t *str)
{
    int n;
    // debug_print("%s", str->s);
    int *s = ksplit(str, '\t', &n);
    if (n != 9) error("Unknown format. %s", str->s);

    char *feature = str->s + s[2];

    int qry = dict_query(G->features, feature);
    if (qry == -1) {
        warnings("Unknown feature type. %s", feature);
        return 1;
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
    int n0, i;
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
            if (!ret) error("Duplicate attribute ? %s", pp->key);
            kh_val((kh_attr_t*)gtf.attr_dict, k) = pair->val == NULL ? NULL : strdup(pair->val);
        }
        free(pp->key);
        if (pp->val) free(pp->val);
    }
    free(pair);
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
    free(s);
    return 0;
}
// key names: gene_id, gene_name, transcript_id,
struct gtf_spec *gtf_read(const char *fname)
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
        if (parse_str(G, &str)) warnings("Skip line %d, %s", line, str.s);
    }
    free(str.s);
    gzclose(fp);
    ks_destroy(ks);
    
    if (G->n_gtf == 0) {
        gtf_destory(G);
        return NULL;
    }

    gtf_build_index(G);
    
    return G;
}
/*
char *gtf_get_gene_name(struct gtf_spec *G, struct gtf_lite *gl)
{
    return G->gene_name->name[gl->gene_name];
}
char *gtf_get_gene_id(struct gtf_spec *G, struct gtf_lite *gl)
{
    return G->gene_id->name[gl->gene_id];
}
char *gtf_get_transcript_id(struct gtf_spec *G, struct gtf_lite *gl)
{
    return G->transcript_id->name[gl->transcript_id];
}
*/
#define idx_start(a) (int)(a>>32)
#define idx_end(a) (int)(a)

static int last_idx = -1;
static int last_id = -1;

void gtf_clean_cache()
{
    last_idx = -1;
    last_id = -1;
}
// if cache == 1, last record will be kept and assume input records have been sorted
struct gtf_lite *gtf_overlap_gene(struct gtf_spec *G, char *name, int start, int end, int *n, int cache)
{
    *n = 0;
    int id = dict_query(G->name, name);
    if (id == -1) return NULL;

    int st = G->ctg[id].idx;
    int ed = st + G->ctg[id].offset-1;
    int ed0 = ed;
    if (end < idx_start(G->idx[st])) return NULL;
    if (start > idx_end(G->idx[ed])) return NULL;

    if (cache == 1 && id == last_id) {
        st = last_idx;
        if (idx_end(G->idx[st]) >= start) goto check_overlap;
    }

    // find the smallest i such that idx_end >= st
    while (st > ed) {
        int mid = st + ((end-start)>>1);
        if (idx_end(G->idx[mid])>=start) st = mid;
        else ed = mid;
    }
    assert(st == ed);

  check_overlap:
    //struct gtf_lite *g0 = &G->gtf[st];
    if (end < G->gtf[st].start) return NULL; // intergenic
    int i;
    int c = 0;
    for (i = st; i <= ed0; ++i) {
        struct gtf_lite *g1 = &G->gtf[i];
        if (g1->start <= end) c++;
        else break;
    }
    *n = c;
    
    last_id = id;
    last_idx = st;
    return &G->gtf[st];
}
void gtf_format_print_test(struct gtf_spec *G)
{
    int i;
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *gl = &G->gtf[i];
        printf("%s\t%d\t%d\t%s\t%s\n", G->name->name[gl->seqname], gl->start, gl->end, G->features->name[gl->type], G->gene_name->name[gl->gene_name]);
        int j;
        for (j = 0; j < gl->n_son; ++j) {
            struct gtf_lite *tx = &gl->son[j];
            printf("%s\t%d\t%d\t%s\t%s\n", G->name->name[tx->seqname], tx->start, tx->end, G->features->name[tx->type], G->gene_name->name[tx->gene_name]);
        }
    }
}

#ifdef GTF_MAIN

int main(int argc, char **argv)
{
    if (argc != 2) error("gtfformat in.gtf");
    struct gtf_spec *G = gtf_read(argv[1]);
    gtf_format_print_test(G);
    gtf_destory(G);
    return 0;
}

#endif
