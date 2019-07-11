#include "utils.h"
#include "bed_lite.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "number.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8193);
KHASH_MAP_INIT_STR(name, int)

static int name2id(kh_name_t *hash, const char *name)
{
    khiter_t k;
    if (hash == NULL) error("Hash is not inited.");    
    k = kh_get(name,hash,name);
    if (k == kh_end(hash)) return -1;
    return kh_val(hash, k);
}
void bed_destroy(struct bedaux *bed)
{
    if (bed == NULL) error("Try to destroy an empty bed.");
    int i;
    for (i = 0; i < bed->n; i++) {
        sfree(bed->names[i]);
        struct bed_chr *c = &bed->c[i];
        int j;
        for (j = 0; j < c->n; ++j) {
            if (c->b[j].name) free(c->b[j].name);
            if (c->b[j].data != NULL) error("Free all extern data first.");
        }
    }
    sfree(bed->c);
    sfree(bed->names);
    kh_destroy(name, bed->name2id);
    sfree(bed);
}
static int parse_str(struct bedaux *bed, kstring_t *str)
{
    int n = 0;
    int *s = ksplit(str, '\t', &n);
    if (s == NULL) return 1;
    if (n < 3) {
        warnings("Bad format.");
        free(s);
        return 1;
    }
    char *chrom = str->s + s[0];
    int tid = name2id((kh_name_t*)bed->name2id, chrom);
    if (tid == -1) {
        if (bed->n == bed->m) {
            bed->m += 10;
            bed->names = realloc(bed->names, sizeof(char*)*bed->m);
            bed->c = realloc(bed->c, sizeof(struct bed_chr)*bed->m);
        }
        bed->names[bed->n] = strdup(chrom);
        struct bed_chr *c = &bed->c[bed->n];
        memset(c, 0, sizeof(struct bed_chr));
        c->id = bed->n;
        tid = bed->n;
        khint_t k;
        int ret;
        k = kh_put(name, (kh_name_t*)bed->name2id, bed->names[bed->n], &ret);
        kh_val((kh_name_t*)bed->name2id, k) = bed->n;
        bed->n++;
    }
    struct bed_chr *c = &bed->c[tid];
    if (c->n == c->m) {
        c->m = c->m == 0 ? 1024 : c->m *2;
        c->b = realloc(c->b, sizeof(struct bed)*c->m);
    }
    struct bed *b = &c->b[c->n];
    memset(b, 0, sizeof(struct bed));
    b->start = str2int(str->s+s[1]);
    b->end = str2int(str->s+s[2]);
    if (b->end < b->start) {
        warnings("The end position located at front of start position.");
        free(s);
        return 1;
    }
    c->n++;
    
    if (n > 3) {
        b->name = strdup(str->s+s[3]);
    }
    else b->name = NULL;

    free(s);
    return 0;
}
struct bedaux *bed_read(const char *fname)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", fname, strerror(errno));
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int line = 0;
    //int m = 100;
    struct bedaux *b = malloc(sizeof(*b));
    memset(b, 0, sizeof(struct bedaux));
    b->name2id = kh_init(name);
    int n = 0;
    while (ks_getuntil(ks, 2, &str, &ret)>=0){
        line++;
        if (str.l == 0) {
            warnings("Line %d is empty. Skip.", line);
            continue;
        }
        if (str.s[0] == '#') continue;
        if (parse_str(b, &str)) {
            str.l = 0;
            continue;
        }
        n++;
    }
    if (str.m) free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    if (n == 0) {
        sfree(b);
        return NULL;
    }
    return b;
}
int bed_select_chrom(struct bedaux *bed, const char *chr)
{
     return name2id((kh_name_t*)bed->name2id, chr);
}
