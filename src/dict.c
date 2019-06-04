#include "utils.h"
#include "dict.h"
#include "htslib/khash.h"

KHASH_MAP_INIT_STR(dict, int)

struct dict *dict_init()
{
    struct dict *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    d->dict = kh_init(dict);
    return d;
}
void dict_destory(struct dict *d)
{
    int i;
    for (i = 0; i < d->n; ++i) free(d->name[i]);
    free(d->name);
    kh_destroy(dict,(kh_dict_t*)d->dict);
    free(d);
}

int dict_query(struct dict *d, char *name)
{
    khint_t k;
    k = kh_get(dict, (kh_dict_t*)d->dict, name);
    if (k == kh_end((kh_dict_t*)d->dict)) return -1;

    return kh_val((kh_dict_t*)d->dict, k);
}

int dict_push(struct dict *d, char *name)
{
    int check = dict_query(d, name);
    if (check == -1) {
        khint_t k;
        int ret;
        k = kh_put(dict, (kh_dict_t*)d->dict, name, &ret);
        if (d->n == d->m) {
            d->m = d->m == 0 ? 32 : d->m<<1;
            d->name = realloc(d->name, sizeof(char *)*d->m);
        }
        d->name[d->n] = strdup(name);
        kh_val((kh_dict_t*)d->dict, k) = d->n;
        check = d->n;
        d->n++;
    }
    return check;
}
