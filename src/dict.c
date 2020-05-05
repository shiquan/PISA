#include "dict.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include <zlib.h>

KHASH_MAP_INIT_STR(name, int)
    
KSTREAM_INIT(gzFile, gzread, 8193);

struct dict {
    int n, m;
    char **name;
    kh_name_t *dict;
    uint32_t *count;
    int assign_value_flag;
    void **value;
};

struct dict *dict_init()
{
    struct dict *D = malloc(sizeof(*D));
    memset(D, 0, sizeof(*D));
    D->dict = kh_init(name);
    return D;
}
void dict_set_value(struct dict *D)
{
    if (D->assign_value_flag == 1) error("Double assign value to a dict.");
    D->assign_value_flag = 1;
    if (D->m > 0) {
        D->value = malloc(D->m*sizeof(void*));
        int i;
        for (i = 0; i < D->m; ++i) D->value[i] = NULL;
    }
}
void *dict_query_value(struct dict *D, int idx)
{
    if (idx < 0 || idx > D->n) return NULL;
    return D->value[idx];
}

void *dict_query_value2(struct dict *D, const char *key)
{
    int idx = dict_query(D, key);
    if (idx == -1) return NULL; // failed to query
    return D->value[idx];
}

static char bit32[4];

void *dict_query_valueInt(struct dict *D, int key)
{
    bit32[0] = (key >> 24) & 0xFF;
    bit32[1] = (key >> 16) & 0xFF;
    bit32[2] = (key >> 8) & 0xFF;
    bit32[3] = key & 0xFF;
    return dict_query_value2(D, bit32);
}

int dict_assign_value(struct dict *D, int idx, void *val)
{
    if (idx < 0 || idx > D->n) return 1;
    D->value[idx] = val;
    return 0;
}

char *dict_name(const struct dict *D, int idx)
{
    assert(idx >= 0 && idx < D->n);
    return D->name[idx];
}
int dict_nameInt(const struct dict *D, int idx)
{
    assert(idx >= 0 && idx < D->n);
    char *name= D->name[idx];
    int x = 0;
    x = x | (name[0] &0xFF);
    x = x << 8;
    x = x | (name[1] &0xFF);
    x = x << 8;
    x = x | (name[2] &0xFF);
    x = x << 8;
    x = x | (name[3] &0xFF);
    return x;
}

int dict_size(const struct dict *D)
{
    return D->n;    
}
uint32_t dict_count(const struct dict *D, int idx)
{
    return D->count[idx];
}
uint32_t dict_count_sum(const struct dict *D)
{
    uint32_t sum = 0;
    int i;
    for (i = 0; i < D->n; ++i) sum+=D->count[i];
    return sum;
}
void dict_destroy(struct dict *D)
{
    int i;
    for (i = 0; i < D->n; ++i) free(D->name[i]);
    // values are actually points, need free pointed values manually
    if (D->assign_value_flag) free(D->value); 
    free(D->name);
    free(D->count);
    kh_destroy(name,D->dict);
    free(D);
}

int dict_query(const struct dict *D, char const *key)
{
    if (key == NULL) error("Trying to query an empty key.");
    khint_t k;
    k = kh_get(name, D->dict, key);
    if (k == kh_end(D->dict)) return -1;
    return kh_val(D->dict, k);
}


int dict_queryInt(const struct dict *D, int key)
{
    bit32[0] = (key >> 24) & 0xFF;
    bit32[1] = (key >> 16) & 0xFF;
    bit32[2] = (key >> 8) & 0xFF;
    bit32[3] = key & 0xFF;

    khint_t k;
    k = kh_get(name, D->dict, bit32);
    if (k == kh_end(D->dict)) return -1;
    return kh_val(D->dict, k);
}

int dict_push0(struct dict *D, char const *key)
{
    if (key == NULL) error("Trying to push an empty key.");
    int ret;
    ret = dict_query(D, key);
    if (ret != -1) {
        D->count[ret]++;
        return ret;
    }
    if (D->n == D->m) {
        D->m = D->m == 0 ? 1024 : D->m<<1;
        D->count = realloc(D->count, sizeof(uint32_t)*D->m);
        assert(D->count);
        int i;
        for (i = D->n; i < D->m; ++i) D->count[i] = 0;
        D->name = realloc(D->name, sizeof(char*)*D->m);
        if (D->assign_value_flag) {
            D->value = realloc(D->value, sizeof(void *)*D->m);
            int j;
            for (j = D->n; j < D->m; ++j) D->value[j] = NULL; // reset  
        }
    }
    D->name[D->n] = strdup(key);
    khint_t k;
    k = kh_put(name, D->dict, D->name[D->n], &ret);
    kh_val(D->dict, k) = D->n;
    return -1;
}

int dict_push(struct dict *D, char const *key)
{
    int ret;
    ret = dict_push0(D, key);
    if (ret != -1) return ret; // already present
    D->count[D->n]++;
    return D->n++;
}
// push new key without increase count
int dict_push1(struct dict *D, char const *key)
{
    int ret;
    ret = dict_push0(D, key);
    if (ret != -1) return ret;
    return D->n++;
}

int dict_pushInt(struct dict *D, int key)
{
    bit32[0] = (key >> 24) & 0xFF;
    bit32[1] = (key >> 16) & 0xFF;
    bit32[2] = (key >> 8) & 0xFF;
    bit32[3] = key & 0xFF;
    return dict_push(D, bit32);
}
int dict_pushInt1(struct dict *D, int key)
{
    bit32[0] = (key >> 24) & 0xFF;
    bit32[1] = (key >> 16) & 0xFF;
    bit32[2] = (key >> 8) & 0xFF;
    bit32[3] = key & 0xFF;
    return dict_push1(D, bit32);
}

int dict_read(struct dict *D, const char *fname)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    while (ks_getuntil(ks, 2, &str, &ret)>=0){
        if (str.l == 0) continue;
        if (str.s[0] == '#') continue;
        if (strcmp(str.s, "Barcode") == 0) {
            warnings("\"Barcode\" in %s looks like a title, skip it. ", fname);
            continue; // emit header
        }
        char *p = str.s;
        char *e = str.s + str.l;
        while (p < e && !isspace(*p)) p++;
        *p = '\0';
        dict_push1(D, str.s); // v0.4, init whitelist but not increase count
    }
    if (str.m) free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return 0;
}

char **dict_names(struct dict *D)
{
    return D->name;
}
// hamming distance
static int check_similar(char *a, char *b, int mis)
{
    int l0, l1;
    l0 = strlen(a);
    l1 = strlen(b);
    assert(l0 == l1);
    int i;
    int m = 0;
    for (i = 0; i < l0; ++i) {
        if (a[i] != b[i]) m++;
        if (m >mis) return 1;
    }
    return 0;
}
// Consider errors during PCR and sequencing, sometime barcodes may contain one or more mismatches.
// this function try to compare each value, allow 1 mismatch, and return the most likely key.
char *dict_most_likely_key(struct dict *D)
{
    uint32_t count = 0;
    char *key = NULL;
    int i;
    for (i = 0; i < D->n; ++i) {
        if (key == NULL) {
            key = D->name[i];
            count = D->count[i];            
        }
        else {
            if (check_similar(key, D->name[i], 1) == 0) {
                if (count < D->count[i]) {
                    count = D->count[i];
                    key = D->name[i];
                }
            }
            else {
                return NULL;
            }
        }            
    }

    return key;
}
