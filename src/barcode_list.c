#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "barcode_list.h"
#include <zlib.h>
KSTREAM_INIT(gzFile, gzread, 8193);
KHASH_MAP_INIT_STR(tag, int)
typedef kh_tag_t taghash_t;

static struct barcode *init_barcode_list(const char *fname, int *n, int *m)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int _m = 100;
    int _n = 0;
    struct barcode *b = malloc(sizeof(struct barcode)*_m);
    while (ks_getuntil(ks, 2, &str, &ret)>=0){
        if (str.l == 0) continue;
        if (str.s[0] == '#') continue;
        if (strcmp(str.s, "Barcode") == 0) continue; // emit header
        char *p = str.s;
        char *e = str.s + str.l;
        while (p < e && !isspace(*p)) p++;
        *p = '\0';
        if (_n == _m) {
            _m = _m<<1;
            b = realloc(b, sizeof(struct barcode)*_m);
        }
        b[_n].s = strdup(str.s);
        b[_n].data = NULL;
        _n++;
    }
    if (str.m) free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    *n = _n;
    *m = _m;
    if (_n == 0) {
        free(b);
        return NULL;
    }
    return b;
}
struct lbarcode *barcode_init()
{
    struct lbarcode *lb = malloc(sizeof(struct lbarcode));
    memset(lb, 0, sizeof(struct lbarcode));
    taghash_t *h = kh_init(tag);
    lb->hash = (void*)h;
    return lb;
}
int barcode_read(struct lbarcode *lb, const char *fname)
{
    lb->b = init_barcode_list(fname, &lb->n, &lb->m);
    khint_t k;
    int i;
    int ret;
    for (i = 0; i < lb->n; ++i) {
        k = kh_get(tag, (taghash_t*)lb->hash, lb->b[i].s);
        if (k == kh_end((taghash_t*)lb->hash)) {
            k = kh_put(tag, (taghash_t*)lb->hash, lb->b[i].s, &ret);
            kh_val((taghash_t*)lb->hash, k) = i;
        }
    }

    return 0;
}
int barcode_select(struct lbarcode *lb, char *s)
{
    taghash_t *hash = (taghash_t*)lb->hash;
    khiter_t k = kh_get(tag, hash, s);
    if (k == kh_end(hash)) return -1;
    return kh_val(hash, k);
}
// for barcode already existed, reture the index
int barcode_push(struct lbarcode *lb, char *s)
{
    int ret = barcode_select(lb, s);
    if (ret != -1) return ret;
    if (lb->n == lb->m) {
        lb->m = lb->m == 0 ? 100 : lb->m<<1;
        lb->b = realloc(lb->b, sizeof(struct barcode)*lb->m);
    }
    struct barcode *b = &lb->b[lb->n];
    b->s =strdup(s);
    b->data = NULL;
    khiter_t k;
    k = kh_put(tag, (taghash_t*)lb->hash, b->s, &ret);
    kh_val((taghash_t*)lb->hash, k) = lb->n;
    lb->n++;
    return lb->n-1;
}
