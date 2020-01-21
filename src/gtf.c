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
    G->ctg = malloc(dict_size(G->name)*sizeof(struct ctg_idx));
    memset(G->ctg, 0, sizeof(struct ctg_idx)*dict_size(G->name));
    G->idx = malloc(G->n_gtf*sizeof(uint64_t));
    int i;
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf_lite *gl = &G->gtf[i];
        G->idx[i] = (uint64_t)gl->start<<32|gl->end;
        G->ctg[gl->seqname].offset++;        
        if (G->ctg[gl->seqname].idx == 0) G->ctg[gl->seqname].idx = i+1;
    }
    for (i = 0; i < dict_size(G->name); ++i) G->ctg[i].idx -= 1; // convert to 0 based
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
    assert(g0->gene_name == gl->gene_name);
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
static int parse_str(struct gtf_spec *G, kstring_t *str, int filter)
{
    int n;
    // debug_print("%s", str->s);
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
        gtf_destory(G);
        return NULL;
    }

    gtf_build_index(G);
    
    return G;
}
void gtf_itr_destory(struct gtf_itr *i)
{
    free(i); 
}

#define idx_start(a) (int)(a>>32)
#define idx_end(a) (int)(a)

struct gtf_itr *gtf_itr_build(struct gtf_spec *G)
{
    struct gtf_itr *i = malloc(sizeof(*i));
    memset(i, 0, sizeof(*i));
    i->G = G;
    i->id = -1;
    return i;
}
// todo: improve performance here
int gtf_query(struct gtf_itr *itr, char *name, int start, int end)
{
    struct gtf_spec *G = itr->G;
    int id = dict_query(G->name, name);
    if (id == -1) return -2; // not this chrom

    int st = G->ctg[id].idx;
    int ed = st + G->ctg[id].offset-1;
    int ed0 = ed;
    if (end < idx_start(G->idx[st])) return -1; // out of range

    // removed because for overlapped transcript, big transcript may cover both ends of short ones
    // if (start > idx_end(G->idx[ed])) return -1; 

    if (id == itr->id) {
        st = itr->st;
        if (idx_end(G->idx[st]) > start) goto check_overlap;
        if (st+1 < ed && idx_start(G->idx[st+1]) > end) return 1; // intergenic
    }

    // FIX: find the smallest i such that start(idx[st]) <= start && start <= end(idx[st])
    // 
    while (st < ed) {
        int mid = st + ((ed-st)>>1);
        if (idx_start(G->idx[mid])>start) st = mid+1;
        else ed = mid;
    }
    if (st != ed) error("%d %d, %d, %d, start : %d, end : %d", st, ed, idx_start(G->idx[st]), idx_end(G->idx[st]), start, end);

  check_overlap:
    if (end < G->gtf[st].start) {
        itr->id = id;
        itr->st = st;
        return 1; // intergenic
    }
    int i;
    int c = 0;
    for (i = st; i <= ed0; ++i) {
        struct gtf_lite *g1 = &G->gtf[i];
        if (g1->start <= end) c++;
        else break;
        // if (c > 4) break; // cover over 4 genes?? impossible
    }

    itr->id = id;
    itr->st = st;
    itr->n = c;
    return 0;
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
