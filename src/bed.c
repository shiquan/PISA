//
#include "utils.h"
#include "dict.h"
#include "region_index.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/vcf.h"
#include "bed.h"
#include "number.h"
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8193)
    
/* struct _ctg_idx { */
/*     int offset; */
/*     int idx; */
/* }; */

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
        if (B->idx)
            region_index_destroy(B->idx[i].idx);
    
    if (B->idx) free(B->idx);
    if (B->ctg) free(B->ctg);
    free(B->bed);
    dict_destroy(B->seqname);
    dict_destroy(B->name);
    free(B);
}
int bed_name2id(struct bed_spec *B, char *name)
{
    return dict_push(B->seqname, name);
}
char* bed_seqname(struct bed_spec *B, int id)
{
    return dict_name(B->seqname, id);
}
void debug_print_bed(struct bed_spec *B)
{
    int i;
    for (i = 0; i < B->n; ++i) {
        fprintf(stderr, "%s\t%d\t%d\n", bed_seqname(B, B->bed[i].seqname), B->bed[i].start, B->bed[i].end);
    }
}
void bed_spec_merge0(struct bed_spec *B, int strand)
{
    qsort(B->bed, B->n, sizeof(struct bed), cmpfunc);
    int i, j;

    for (i = 0; i < B->n; ++i) {
        struct bed *bed0 = &B->bed[i];
        for (j = i+1; j < B->n; ++j) {
            struct bed *bed = &B->bed[j];
            if (bed->seqname == -1) continue;
            if (bed0->seqname == -1) { // last record has been reset
                memcpy(bed0, bed, sizeof(struct bed));
                bed->seqname = -1; //reset
                continue;
            }
            if (bed->seqname != bed0->seqname) {
                memmove(B->bed+i+1, B->bed+j, (B->n-j)*sizeof(struct bed));
                break;
            }
            // check strand sensitive
            if (strand == 1 && bed->strand != bed0->strand) continue;
            
            if (bed->start >= bed0->end) break; // nonoverlap, move to next record
            if (bed0->end < bed->end) bed0->end = bed->end; // enlarge overlapped region
            
            bed->seqname = -1; // reset this record
        }

        if (bed0->seqname == -1) break; // last record
    }

    for (i = B->n-1; i >= 0; --i) {
        struct bed *bed = &B->bed[i];
        if(bed->seqname!= -1) break;
    }
    
    B->n = i+1;
}
static void bed_build_index(struct bed_spec *B)
{
    qsort(B->bed, B->n, sizeof(struct bed), cmpfunc);

    B->ctg = malloc(dict_size(B->seqname)*sizeof(struct _ctg_idx));
    //memset(B->ctg, 0, sizeof(struct _ctg_idx)*dict_size(B->seqname));    
    B->idx = malloc(dict_size(B->seqname)*sizeof(struct bed_idx));
    
    int i;
    for (i = 0; i < dict_size(B->seqname); ++i) {
        // init ctg
        B->ctg[i].offset = 0;
        B->ctg[i].idx = -1;
        // init idx
        B->idx[i].idx = region_index_create();
    }

    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        B->ctg[bed->seqname].offset++;
        if (B->ctg[bed->seqname].idx == -1) B->ctg[bed->seqname].idx = i; // 0-based
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

    if (n >= 4) {
        char *name = str->s+s[3];
        if (strlen(name) == 1) {
            if (name[0] == '.' || name[0] == '*') bed->name = -1;
            else  bed->name = dict_push(B->name, name);
        } else {
            bed->name = dict_push(B->name, name);  
        } 
    }

    bed->strand = -1; // undefined
    if (n >= 6) {
        char *strand = str->s+s[5];
        if (*strand == '+') bed->strand = BED_STRAND_FWD;
        if (*strand == '-') bed->strand = BED_STRAND_REV;
    }
    B->n++;

    free(s);
    return 0;
}

int bed_spec_push(struct bed_spec *B, struct bed *bed)
{
    if (B->n == B->m) {
        B->m = B->m == 0 ? 32 : B->m<<1;
        B->bed = realloc(B->bed, sizeof(struct bed)*B->m);
    }
    struct bed *bed0 = &B->bed[B->n];
    memcpy(bed0, bed, sizeof(struct bed));
    return B->n++;
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
static struct var *var_init()
{
    struct var *v = malloc(sizeof(*v));
    v->ref = malloc(sizeof(kstring_t));
    v->ref->l = v->ref->m = 0;
    v->ref->s = NULL;
    v->alt = malloc(sizeof(kstring_t));
    v->alt->l = v->alt->m = 0;
    v->alt->s = NULL;
    return v;
}

struct bed_spec *bed_read_vcf(const char *fn)
{
    htsFile *fp = hts_open(fn, "r");
    if (fp == NULL) error("%s : %s.", fn, strerror(errno));

    htsFormat type = *hts_get_format(fp);
    if (type.format != vcf && type.format != bcf)
        error("Unsupport input format. -vcf only accept VCF/BCF file.");

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (hdr == NULL) error("Failed parse header of input.");

    struct bed_spec *B = bed_spec_init();
    
    bcf1_t *v = bcf_init();
    
    while(bcf_read(fp, hdr, v) == 0) {
        if (v->rid == -1) continue;
        bcf_unpack(v, BCF_UN_INFO);
        const char *name = bcf_hdr_int2id(hdr, BCF_DT_CTG, v->rid);
        int seqname = dict_push(B->seqname, name);
        if (B->n == B->m) {
            B->m = B->m == 0 ? 32 : B->m<<1;
            B->bed = realloc(B->bed, sizeof(struct bed)*B->m);
        }
        struct bed *bed = &B->bed[B->n];
        bed->seqname = seqname;
        bed->start = v->pos;
        bed->end = v->pos+v->rlen;
        bed->name = -1;
        bed->strand = -1;
        struct var *var = var_init();
        kputs(v->d.allele[0], var->ref);
        if (v->n_allele > 1) kputs(v->d.allele[1], var->alt);

        bed->data = var;
        B->n++;
    }
    bcf_destroy(v);
    hts_close(fp);
    bcf_hdr_destroy(hdr);

    bed_build_index(B);
    
    return B;
}

void bed_spec_var_destroy(struct bed_spec *B)
{
    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        struct var *var = (struct var*)bed->data;
        if (var->ref->m) free(var->ref->s);
        free(var->ref);
        if (var->alt->m) free(var->alt->s);
        free(var->alt);
        free(var);
    }
    bed_spec_destroy(B);
}

struct region_itr *bed_query(const struct bed_spec *B, char *name, int start, int end, int strand)
{
    int id = dict_query(B->seqname, name);
    if (id == -1) return NULL;

    if (start < 0) start = 0;
    if (end < start) {
        warnings("Bad ranger, %s:%d-%d", name, start, end);
        return NULL;
    }
    
    int st = B->ctg[id].idx; // 0 based
    if (end < B->bed[st].start) return NULL; // out of range

    struct region_index *idx = B->idx[id].idx;
    struct region_itr *itr = region_query(idx, start, end);
    if (itr == NULL) return NULL;
    int i;
    for (i = 0; i < itr->n;) {
        struct bed *bed = itr->rets[i];
        if (bed->start >= end || bed->end < start) { // start is 0 base
            memmove(itr->rets+i, itr->rets+i+1, (itr->n-i-1)*sizeof(void*));
            itr->n--;
            continue;
        }
        if (strand != BED_STRAND_IGN) { // check strand
            if (strand != bed->strand) {
                memmove(itr->rets+i, itr->rets+i+1, (itr->n-i-1)*sizeof(void*));
                itr->n--;
                continue;
            }
        }
        i++;
    }
    if (itr->n == 0) {
        region_itr_destroy(itr);
        return NULL;
    }
    qsort((struct bed**)itr->rets, itr->n, sizeof(struct bed*), cmpfunc1);

    return itr;
}
// return 0 on nonoverlap, 1 on overlap
int bed_check_overlap(const struct bed_spec *B, char *name, int start, int end, int strand)
{
    struct region_itr *itr = bed_query(B, name, start, end, strand);
    if (itr == NULL) return 0;
    region_itr_destroy(itr);
    return 1;
}
