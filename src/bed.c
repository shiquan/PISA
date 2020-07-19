//
#include "utils.h"
#include "dict.h"
#include "region_index.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "bed.h"
#include "number.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8193)
    
struct _ctg_idx {
    int offset;
    int idx;
};

struct bed_idx {
    struct region_index *idx;
};

static int cmpfunc(const void *_a, const void *_b)
{
    struct bed *a = (struct bed*)_a;
    struct bed *b = (struct bed*)_b;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    return a->end - b->end;
}
static int cmpfunc1(const void *_a, const void *_b)
{
    struct bed *a = *(struct bed **)_a;
    struct bed *b = *(struct bed **)_b;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    return a->end - b->end;    
}

struct bed_spec *bed_spec_init()
{
    struct bed_spec *B = malloc(sizeof(*B));
    memset(B, 0, sizeof(*B));
    B->seqname = dict_init();
    B->name    = dict_init();
    return B;
}

void bed_spec_destroy(struct bed_spec *B)
{
    int i;
    for (i = 0; i < dict_size(B->seqname); ++i)
        region_index_destroy(B->idx[i].idx);
    free(B->idx);
    free(B->ctg);
    free(B->bed);
    dict_destroy(B->seqname);
    dict_destroy(B->name);
}

static void bed_build_index(struct bed_spec *B)
{
    qsort(B->bed, B->n, sizeof(struct bed), cmpfunc);

    B->ctg = malloc(dict_size(B->seqname)*sizeof(struct _ctg_idx));
    memset(B->ctg, 0, sizeof(struct _ctg_idx)*dict_size(B->seqname));
    
    B->idx = malloc(dict_size(B->seqname)*sizeof(struct bed_idx));
    
    int i;
    for (i = 0; i < dict_size(B->seqname); ++i)
        B->idx[i].idx = region_index_create();

    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        B->ctg[bed->seqname].offset++;
        if (B->ctg[bed->seqname].idx == 0) B->ctg[bed->seqname].idx = i+1;
        index_bin_push(B->idx[bed->seqname].idx, bed->start, bed->end, bed);
    }
}

static int parse_str(struct bed_spec *B, kstring_t *str)
{
    int n;
    int *s = ksplit(str, '\t', &n);
    if (n < 3) error("Unknown format. %s", str->s);

    int seqname = dict_push(B->seqname, str->s + s[0]);
    int start   = str2int(str->s+s[1]);
    int end     = str2int(str->s+s[2]);
    if (start == 0 && end == 0) return 1;

    if (end < start) error("Bad range: %s:%d-%d", str->s+s[0], start, end);

    if (B->n == B->m) {
        B->m = B->m == 0 ? 32 : B->m<<1;
        B->bed = realloc(B->bed, sizeof(struct bed)*B->m);
    }
    struct bed *bed = &B->bed[B->n];
    bed->seqname = seqname;
    bed->start   = start;
    bed->end     = end;
    bed->name    = -1;
    if (n >= 4) bed->name = dict_push(B->name, str->s+s[3]);

    bed->strand = -1; // undefined
    if (n >= 6) {
        char *strand = str->s+s[5];
        if (*strand == '+') bed->strand = 0;
        if (*strand != '-') bed->strand = 1;
    }
    B->n++;

    free(s);
    return 0;
}

struct bed_spec *bed_read(const char *fname)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    if (fp == NULL) error("%s : %s.", fname, strerror(errno));

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    struct bed_spec *B = bed_spec_init();
    int line = 0;
    
    while (ks_getuntil(ks, 2, &str, &ret)>=0) {
        line ++;
        if (str.l == 0) {
            warnings("Line %d is empty. Skip", line);
            continue;
        }
        if (str.s[0] == '#') continue;
        parse_str(B, &str);
    }

    free(str.s);
    gzclose(fp);
    ks_destroy(ks);

    if (B->n == 0) {
        bed_spec_destroy(B);
        return NULL;
    }

    bed_build_index(B);
    return B;
}

struct region_itr *bed_query(struct bed_spec *B, char *name, int start, int end)
{
    int id = dict_query(B->seqname, name);
    if (id == -1) return NULL;

    if (start < 0) start = 0;
    if (end < start) {
        warnings("Bad ranger, %s:%d-%d", name, start, end);
        return NULL;
    }

    int st = B->ctg[id].idx-1; // 0 based
    if (end < B->bed[st].start) return NULL; // out of range

    struct region_index *idx = B->idx[id].idx;
    struct region_itr *itr = region_query(idx, start, end);
    if (itr == NULL) return NULL;
    int i;
    for (i = 0; i < itr->n;) {
        struct bed *bed = itr->rets[i];
        if (bed->start > end || bed->end < start) {
            memmove(itr->rets+i, itr->rets+i+1, itr->n-i-1);
            itr->n--;
        }
        else i++;
    }
    if (itr->n == 0) {
        region_itr_destroy(itr);
        return NULL;
    }
    qsort((struct bed**)itr->rets, itr->n, sizeof(struct bed*), cmpfunc1);

    return itr;
}
// return 0 on nonoverlap, 1 on overlap
int bed_check_overlap(struct bed_spec *B, char *name, int start, int end)
{
    struct region_itr *itr = bed_query(B, name, start, end);
    if (itr == NULL) return 0;
    region_itr_destroy(itr);
    return 1;
}
