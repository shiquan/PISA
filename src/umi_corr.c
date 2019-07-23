// Correct UMI for each cell barcode at one gene

#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "number.h"
#include "umi_corr.h"

/*

  UMI structure

 */
struct kmer {
    int cnt;
    int idx; // correct to idx
    char *kmer;
};

KHASH_MAP_INIT_STR(name, int)

static int check_similar(char *a, char *b)
{
    int l1, l2;
    l1 = strlen(a);
    l2 = strlen(b);
    if (l1 == 0 || l2 == 0) error("Try to compare an empty string");
    if (l1 != l2) error("Try to compare two unequal string.");
    int i, m = 0;
    for (i = 0; i < l1; ++i) {
        if (a[i] != b[i]) m++;
        if (m > 1) break;
    }
    if (m > 1) return 1;
    return 0;
}

struct umi_tag {
    int n, m;
    int l;
    struct kmer *kmers;
    void *dict;
};

struct umi_tag *umi_tag_build()
{
    struct umi_tag *U = malloc(sizeof(*U));
    memset(U, 0, sizeof(*U));
    U->dict = kh_init(name);
    return U;
}

void umi_tag_push(struct umi_tag *U, char const *_s)
{
    khint_t k;
    int l;
    char *s = strdup(s);
    l = strlen(s);
    if (l == 0) error("Try to push empty string.");
    if (U->l == 0) U->l = l;
    else if (U->l != l) error("Try to push an unequal kmers.");
    int ret;
    k = kh_put(name, (kh_name_t*)U->dict, s, &ret);
    if (ret) {
        if (U->n == U->m) {
            U->m = U->m+10;
            U->kmers = realloc(U->kmers, U->m*sizeof(struct kmer));
        }
        struct kmer *mer = &U->kmers[U->n];
        kh_val((kh_name_t*)U->dict, k) = U->n;
        mer->cnt = 1;
        mer->idx = U->n;
        mer->kmer = s;
        U->n++;
    }
    else {
        int idx = kh_val((kh_name_t*)U->dict, k);
        U->kmers[idx].cnt++;
        free(s);
    }
}
void umi_tag_corr(struct umi_tag *U)
{
    if (U->n == 1) return;
    int i;
    for (i = 0; i < U->n; ++i) {
        int j;
        for (j = i+1; j < U->n; ++j) {
            struct kmer *k1 = &U->kmers[i];
            struct kmer *k2 = &U->kmers[j];
            if (check_similar(k1->kmer, k2->kmer) == 0) {
                if (k1->cnt < k2->cnt) k1->idx = j;
                else k2->idx = i;
            }
        }
    }
    // redirect the idx    
    for (i = 0; i < U->n; ++i) {
        struct kmer *k1 = &U->kmers[i];
        int idx = k1->idx;
        do {
            struct kmer *k2 = &U->kmers[idx];
            if (k2->idx != idx) idx = k2->idx;
            else break;
        } while (1);
        k1->idx = idx;
    }    
}
char *umi_tag_query(struct umi_tag *U, char const *s)
{
    assert(s != NULL);
    int l;
    l = strlen(s);
    if (l != U->l) error("Trying to query a sequence with different length.");
    khint_t k;
    k = kh_get(name, (kh_name_t*)U->dict, s);
    if (k == kh_end((kh_name_t*)U->dict)) return NULL;
    int idx = kh_val((kh_name_t*)U->dict, k);
    return U->kmers[idx].kmer;
}

void umi_tag_clear(struct umi_tag *U)
{
    int i;
    for (i = 0; i < U->n; ++i)
        free(U->kmers[i].kmer);
    free(U->kmers);
    kh_destroy(name, (kh_name_t*)U->dict);
    memset(U, 0, sizeof(*U));
}
void umi_tag_destory(struct umi_tag *U)
{
    umi_tag_clear(U);
    free(U);
}

/*
  Tag correct structure
 */

struct block_umi {
    char *name;
    struct umi_tag *umi;
};

/*
struct corr_tag {
    int n, m;
    struct block_umi *umi;
    void *dict;
};
*/

struct corr_tag *corr_tag_build()
{
    struct corr_tag *C = malloc(sizeof(*C));
    memset(C, 0, sizeof(*C));
    C->dict = kh_init(name);
    return C;
}
void corr_tag_destory(struct corr_tag *C)
{
    int i;
    for (i = 0; i < C->n; ++i) {
        free(C->umi[i].name);
        umi_tag_destory(C->umi[i].umi);
    }
    kh_destroy(name, (kh_name_t*)C->dict);
}
void corr_tag_push(struct corr_tag *C, char const *n, char const *u)
{
    khint_t k;
    k = kh_get(name, (kh_name_t*)C->dict, n);
    if (k == kh_end((kh_name_t*)C->dict)) {
        if (C->n == C->m) {
            C->m = C->m == 0 ? 1024 : C->m*2;
            C->umi = realloc(C->umi, C->m*sizeof(struct block_umi));
        }
        struct block_umi *umi = &C->umi[C->n++];
        umi->name = strdup(n);
        umi->umi = umi_tag_build();
        umi_tag_push(umi->umi, u);
    }
    else {
        int idx = kh_val((kh_name_t*)C->dict, k);
        struct block_umi *umi = &C->umi[idx];
        umi_tag_push(umi->umi, u);
    }
}

char *corr_tag_retrieve(struct corr_tag *C, char const *n, char const *u)
{
    khint_t k;
    k = kh_get(name, (kh_name_t*)C->dict, n);
    if (k == kh_end((kh_name_t*)C->dict)) return NULL;
    int idx = kh_val((kh_name_t*)C->dict, k);
    return umi_tag_query(C->umi[idx].umi, u);
}
