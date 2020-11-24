//
#include "utils.h"
#include "htslib/khash.h"
#include "region_index.h"

struct binlist {
    int n, m;
    void **a;
};

KHASH_MAP_INIT_INT(bin, struct binlist)

struct region_index {
    khash_t(bin) *idx;
    int size;
};

struct region_index *region_index_create()
{
    struct region_index *idx = malloc(sizeof(struct region_index));
    idx->idx = kh_init(bin);
    idx->size = 0;
    return idx;
}

void region_index_destroy(struct region_index *idx)
{
    khint_t k;
    for (k = kh_begin(idx->idx); k != kh_end(idx->idx); ++k) {
        if (kh_exist(idx->idx, k)) {
            free(kh_val(idx->idx, k).a);
        }
    }
    kh_destroy(bin,idx->idx);
    free(idx);
    idx=NULL;
}

// copied from tabix/index.c
static inline int ti_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

void index_bin_push(struct region_index *idx, uint32_t start, uint32_t end, void *new)
{
    khint_t k;
    int ret;
    struct binlist *l;
    int bin = ti_reg2bin(start, end);
    k = kh_put(bin, idx->idx, bin, &ret);
    l = &kh_val(idx->idx, k);
    if (ret) { // not present
        l->m = 1; l->n = 0;
        l->a = calloc(1, sizeof(void*));
        idx->size++;
    }
    if (l->m == l->n) {
        l->m <<=1;
        l->a = realloc(l->a, l->m*sizeof(void*));
    }
    l->a[l->n++] = new;
}

#define MAX_BIN 37450 // =(8^6-1)/7+1

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
{
	int i = 0, k;
	if (beg >= end) return 0;
	if (end >= 1u<<29) end = 1u<<29;
	--end;
	list[i++] = 0;
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
}
// start is 0 based, end is 1 based
struct region_itr *region_query(struct region_index *idx, int start, int end)
{
    if (start < 0) start = 0;
    if (end < start) return NULL;
    if (idx->size == 0) return NULL;
    int n=0, i, n_bin;
    khint_t k;
    uint16_t *bins  = (uint16_t*)calloc(MAX_BIN, 2);
    n_bin = reg2bins(start, end, bins);

    for (i = 0; i < n_bin; ++i) 
        if ((k = kh_get(bin, idx->idx, bins[i])) != kh_end(idx->idx))
            n += kh_val(idx->idx, k).n;
    
    if (n == 0) {
        free(bins);
        return NULL;
    }

    struct region_itr *itr = malloc(sizeof(*itr));
    itr->n = n;
    itr->rets = malloc(sizeof(void*)*n);
    int l = 0;
    for (i = 0; i < n_bin; ++i)
        if ((k = kh_get(bin, idx->idx, bins[i])) != kh_end(idx->idx)) {
            int j;
            struct binlist *list = &kh_val(idx->idx, k);
            for (j = 0; j < list->n; ++j)
                itr->rets[l++] = list->a[j];                
        }

    free(bins);

    return itr;
}

void region_itr_destroy(struct region_itr *itr)
{
    free(itr->rets);
    free(itr);
    itr=NULL;
}
