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

static const char *bed_anno_type[] = {
    "unknown",
    "promoter", // use for old version, now be replaced by upstream
    "utr3",
    "utr5",
    "exon",
    "multiexons",
    "exonintron",
    "whole_gene",
    "multigenes",
    "intron",
    "antisense_utr3",
    "antisense_utr5",
    "antisense_exon",
    "antisense_intron",
    "antisense_complex",
    "flank",  // flank region
    "upstream",
    "downstream",
    "antisense_upstream",
    "antisense_downstream",
    "intergenic",
    "unknown_chromosome",
};

const char *bed_typename(int type)
{
    return bed_anno_type[type];
}

struct bed_ext *bed_ext_init()
{
    struct bed_ext *e = malloc(sizeof(*e));
    memset(e, 0, sizeof(*e));
    return e;
}

struct bed_idx {
    struct region_index *idx;
};

void bed_copy(struct bed *a, struct bed *b)
{
    a->seqname = b->seqname;
    a->name = b->name;
    a->start = b->start;
    a->end = b->end;
    a->strand = b->strand;
    a->data = b->data;
}
struct bed_spec *bed_spec_dup(struct bed_spec *B0)
{
    struct bed_spec *B = bed_spec_init();
    B->seqname = dict_dup(B0->seqname);
    B->name = dict_dup(B0->name);
    B->idx = NULL;
    B->n = B->m = B0->n;
    B->bed = malloc(B->n*sizeof(struct bed));
    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        struct bed *bed1 = &B0->bed[i];
        bed_copy(bed, bed1);
    }
    // B->ext = B0->ext;
    return B;
}
static int cmpfunc(const void *_a, const void *_b)
{
    struct bed *a = (struct bed*)_a;
    struct bed *b = (struct bed*)_b;
    if (a->seqname == -1 && b->seqname != -1) return 1;
    if (b->seqname == -1 && a->seqname != -1) return -1;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    if (a->end != b->end) return a->end - b->end;
    if (a->strand != b->strand) return a->strand - b->strand;
    return a->name - b->name;
}
static int cmpfunc1(const void *_a, const void *_b)
{
    struct bed *a = *(struct bed **)_a;
    struct bed *b = *(struct bed **)_b;
    if (a->seqname != b->seqname) return a->seqname - b->seqname;
    if (a->start != b->start) return a->start - b->start;
    if (a->end != b->end) return a->end - b->end;
    if (a->strand != b->strand) return a->strand - b->strand;
    return a->name - b->name;
}

struct bed_spec *bed_spec_init()
{
    struct bed_spec *B = malloc(sizeof(*B));
    memset(B, 0, sizeof(*B));
    B->seqname = dict_init();
    B->name    = dict_init();
    return B;
}

void bed_ext_destroy(struct bed_ext *e)
{
    if (e) {
        if (e->genes) {
            int i;
            for (i = 0; i < e->n; ++i)
                if (e->genes[i]) free(e->genes[i]);
            free(e->genes);
        }
        free(e);
    }
}
void bed_spec_ext_destroy(struct bed_spec *B)
{
    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed_ext *e = (struct bed_ext *)B->bed[i].data;
        if (e) bed_ext_destroy(e);
    }
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

void bed_spec_shrink(struct bed_spec *B, int up, int down)
{
    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        if (strand_is_minus(bed->strand)) {
            bed->start = bed->start +down;
            bed->end = bed->end - up;
        } else {
            bed->start = bed->start + up;
            bed->end = bed->end - down;
        }
        if (bed->start < 0) bed->start = 0;
        if (bed->end <= bed->start) {
            bed->start = 0;
            bed->end = 0;
        }
    }
}
void bed_spec_expand(struct bed_spec *B, int up, int down)
{
    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        if (strand_is_minus(bed->strand)) {
            bed->start = bed->start - down;
            if (bed->start < 0) bed->start = 0;
            bed->end = bed->end + up;
        } else {
            bed->start = bed->start - up;
            if (bed->start < 0) bed->start = 0;
            bed->end = bed->end + down;
        }
    }
}
void bed_spec_remove0(struct bed_spec *B)
{
    int i;
    for (i = 0; i < B->n;) {
        struct bed *bed = &B->bed[i];
        if (bed->start == 0 && bed->end == 0) {
            if (B->n == i+1) {
                B->n = 0; // empty
                break;
            }
            memmove(B->bed+i, B->bed+i+1, (B->n-i-1)*sizeof(struct bed));
            B->n--;
            continue;
        }
        ++i;
    }
}
void bed_spec_sort(struct bed_spec *B)
{
    qsort(B->bed, B->n, sizeof(struct bed), cmpfunc);    
}
void bed_spec_dedup(struct bed_spec *B, int check_name)
{
    bed_spec_sort(B);

    int i;
    for (i = 0; i < B->n; i++) {
        struct bed *bed0 = &B->bed[i];
        if (bed0->seqname == -1) continue;
        for (int j = i + 1; j < B->n; ++j) {
            struct bed *bed1 = &B->bed[j];
            if (bed1->seqname == -1) continue;
            if (check_name == 1) {
                if (bed0->seqname == bed1->seqname && bed0->start == bed1->start &&
                    bed0->end == bed1->end &&
                    bed0->strand == bed1->strand && bed0->name == bed1->name) {
                    bed1->seqname = -1;
                    continue;
                }
            } else {
                if (bed0->seqname == bed1->seqname && bed0->start == bed1->start &&
                    bed0->end == bed1->end &&
                    bed0->strand == bed1->strand) {
                    bed1->seqname = -1;
                    continue;
                }
            }

            if (bed1->seqname != bed0->seqname || bed1->start != bed0->start || bed1->end != bed0->end)
                break;
        }
    }
    bed_spec_sort(B);
    for (i = B->n-1; i >= 0; --i) {
        struct bed *bed = &B->bed[i];
        if (bed->seqname != -1) break;
    }

    B->n = ++i;
}
void bed_spec_merge0(struct bed_spec *B, int strand, int check_name)
{
    bed_spec_sort(B);
    int i, j;

    for (i = 0; i < B->n; ++i) {
        struct bed *bed0 = &B->bed[i];
        if (bed0->seqname == -1) continue;
        
        if (strand == 0) bed0->strand = BED_STRAND_UNK; // reset strand

        for (j = i+1; j < B->n; ++j) {
            struct bed *bed = &B->bed[j];
            if (bed->seqname == -1) continue;

            if (bed->seqname != bed0->seqname) {
                break;
            }
            // check strand sensitive
            if (strand == 1 && bed->strand != bed0->strand) continue;
            if (check_name && bed->name != -1 && bed0->name != -1 && bed->name != bed0->name) continue;
            if (bed->start > bed0->end) break; // nonoverlap, move to next record
            if (bed0->end < bed->end) bed0->end = bed->end; // enlarge overlapped region
            
            bed->seqname = -1; // reset this record
        }
    }

    bed_spec_sort(B);
    
    for (i = B->n-1; i >= 0; --i) {
        struct bed *bed = &B->bed[i];
        if(bed->seqname!= -1) break;
    }
    
    B->n = i+1;
}
void bed_spec_merge1(struct bed_spec *B, int strand, int up, int down, int min_length, int check_name)
{
    bed_spec_expand(B, up, down);
    bed_spec_merge0(B, strand, check_name);    
}
void bed_spec_merge2(struct bed_spec *B, int strand, int gap, int min_length, int check_name)
{
    bed_spec_merge0(B, strand, check_name);
    if (gap > 0) {
        bed_spec_expand(B, 0, gap);
        bed_spec_shrink(B, 0, gap);
    }

    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        if (bed->end - bed->start < min_length) {
            bed->seqname = -1;
        }
    }
    bed_spec_merge0(B, strand, check_name);
}

void bed_build_index(struct bed_spec *B)
{
    qsort(B->bed, B->n, sizeof(struct bed), cmpfunc);

    B->ctg = malloc(dict_size(B->seqname)*sizeof(struct _ctg_idx));
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

int bed_spec_push0(struct bed_spec *B, const char *seqname, int start, int end, int strand, const char *name, void *ext)
{
    if (seqname == NULL) error("No seqname.");
    if (start < 0 || end < 0) error("Unset start or end position.");
    if (strand != 1 && strand != 0 && strand != -1) error("Unknown strand.");
    
    if (B->n == B->m) {
        B->m = B->m == 0 ? 32 : B->m<<1;
        B->bed = realloc(B->bed, sizeof(struct bed)*B->m);
    }
    
    int seqname_id = dict_push(B->seqname, seqname);
    int name_id = -1;
    if (name) name_id = dict_push(B->name, name);
    
    struct bed *b = &B->bed[B->n];
    b->seqname = seqname_id;
    b->start = start;
    b->end = end;
    b->strand = strand;
    b->name = name_id;
    b->data = ext;
    return B->n++;
}

int bed_spec_push1(struct bed_spec *B, int seqname, int start, int end, int strand, int name, void *ext)
{
    if (B->n == B->m) {
        B->m = B->m == 0 ? 32 : B->m<<1;
        B->bed = realloc(B->bed, sizeof(struct bed)*B->m);
    }
    
    struct bed *b = &B->bed[B->n];
    b->seqname = seqname;
    b->start = start;
    b->end = end;
    b->strand = strand;
    b->name = name;
    b->data = ext;
    return B->n++;
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

struct bed_spec *bed_read0(struct bed_spec *B, const char *fname)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    if (fp == NULL) error("%s : %s.", fname, strerror(errno));

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
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
    return B;
}

struct bed_spec *bed_read(const char *fname)
{
    struct bed_spec *B = bed_spec_init();

    B = bed_read0(B, fname);

    if (B->n == 0) {
        bed_spec_destroy(B);
        return NULL;
    }
    
    bed_build_index(B);
    return B;
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
    B->ext = hdr;
    
    bcf1_t *v = bcf_init();
    
    while (bcf_read(fp, hdr, v) == 0) {
        if (v->rid == -1) continue;
        bcf_unpack(v, BCF_UN_STR);
        if (v->n_allele == 1) continue;
        if (v->n_allele == 2) {
            if (v->d.allele[1][0] == '<') continue; // skip <NON_REF>, <X>, <*>
        }
        
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
        bcf1_t *var = bcf_init();
        bed->data = bcf_copy(var, v);
        bcf_unpack(var, BCF_UN_STR);
        bcf_unpack(var, BCF_UN_INFO);
        B->n++;
    }
    bcf_destroy(v);
    hts_close(fp);

    bed_build_index(B);
    return B;
}

void bed_spec_var_destroy(struct bed_spec *B)
{
    bcf_hdr_t *hdr = (bcf_hdr_t*)B->ext;
    if (hdr == NULL) error("No vcf hdr.");
    bcf_hdr_destroy(hdr);

    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];
        bcf1_t *v = (bcf1_t*)bed->data;
        bcf_destroy(v);
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
// chrom, start, end, name, score[reserved], strand,
// ext: n_gene, gene(s), functional type, nearby gene for integenic, nearby distance
void bed_spec_write0(struct bed_spec *B, FILE *out, int ext, int gene_as_name)
{
    int i;
    for (i = 0; i < B->n; ++i) {
        struct bed *bed = &B->bed[i];

        if (bed->seqname == -1) continue;
        
        if (ext) {
            struct bed_ext *e = (struct bed_ext*)bed->data;
            if (e) {
                kstring_t name = {0,0,0};
                if (e->n == 1) {
                    kputs(e->genes[0], &name);  
                } else {
                    if (e->n > 3) {
                        ksprintf(&name, "%s,...,%s", e->genes[0], e->genes[e->n-1]);
                    } else {
                        int k;
                        for (k = 0; k < e->n; ++k) {
                            if (k) kputc(',', &name);
                            kputs(e->genes[k], &name);
                        }
                    }
                }
                
                if (gene_as_name && name.m > 0) {
                    fprintf(out, "%s\t%d\t%d\t%s\t.\t%c", dict_name(B->seqname, bed->seqname),
                            bed->start, bed->end, name.s, //bed->name == -1 ? "." : dict_name(B->name, bed->name),
                            ".+-"[bed->strand+1]
                        );
                } else {
                    fprintf(out, "%s\t%d\t%d\t%s\t.\t%c", dict_name(B->seqname, bed->seqname),
                            bed->start, bed->end, bed->name == -1 ? "." : dict_name(B->name, bed->name),
                            ".+-"[bed->strand+1]
                        );                
                }
                
                fprintf(out, "\t%d\t", e->n);

                if (name.l)fputs(name.s, out);
                else fputc('.', out);
                
                fprintf(out, "\t%s", bed_typename(e->type));
                
                if (name.m) free(name.s);
            } else {
                fprintf(out, "%s\t%d\t%d\t%s\t.\t%c\t0\t.\tintergenic", dict_name(B->seqname, bed->seqname),
                        bed->start, bed->end, bed->name == -1 ? "." : dict_name(B->name, bed->name),
                        ".+-"[bed->strand+1]
                    );
            }
        } else {
            fprintf(out, "%s\t%d\t%d\t%s\t.\t%c", dict_name(B->seqname, bed->seqname),
                    bed->start, bed->end, bed->name == -1 ? "." : dict_name(B->name, bed->name),
                    ".+-"[bed->strand+1]
                );
        }
        
        fputc('\n', out);
    }
}

void bed_spec_write(struct bed_spec *B, const char *fn, int ext, int gene_as_name)
{
    FILE *fp;
    if (fn == NULL) fp = stdout;
    else {
        fp = fopen(fn, "w");
        if (fp == NULL) error("%s : %s", fn, strerror(errno));
    }

    bed_spec_write0(B, fp, ext, gene_as_name);

    fclose(fp);
}
void bed_spec_seqname_from_bam(struct bed_spec *B, bam_hdr_t *hdr)
{
    if (dict_size(B->seqname) > 0) error("Init non-empty seqnames.");
    int i;
    for (i = 0; i < hdr->n_targets; ++i)
        dict_push(B->seqname, hdr->target_name[i]);
}
void bed_spec_seqname_from_gtf(struct bed_spec *B, struct gtf_spec *G)
{
    if (dict_size(B->seqname) > 0) error("Init non-empty seqnames.");
    int i;
    for (i = 0; i < dict_size(G->name); ++i)
        dict_push(B->seqname, dict_name(G->name, i));
}
static int cmpint(const void *_a, const void *_b)
{
    return (*(int*)_a - *(int*)_b);
}
struct bed_spec *bed_spec_flatten(struct bed_spec *B0, int offset)
{
    struct bed_spec *B = bed_spec_dup(B0);
    bed_spec_dedup(B, 1);

    struct bed_spec *new = bed_spec_dup(B);
    new->n = 0;

    if (offset != 0) offset = 1;
    
    int n = 0, m = 0;
    int *bps = NULL;
    
    for (int i = 0; i < B->n; ++i) {
        n = 0;
        struct bed *bed0 = &B->bed[i];
        if (bed0->seqname == -1) continue;
        for (int j = i+1; j < B->n; ++j) {
            struct bed *bed = &B->bed[j];
            if (bed->seqname == -1) continue;
            if (bed->seqname != bed0->seqname) break;

            if (bed->strand != bed0->strand) continue;
            if (bed->name != bed0->name) continue;
            if (bed->start > bed0->end) break;
            if (bed->start > bed0->start) {
                // push start
                if (n == m) {
                    m = m == 0 ? 4 : m *2;
                    bps = realloc(bps, sizeof(int)*m);
                }
                bps[n++] = bed->start;
            }

            if (n == m) {
                m = m == 0 ? 4 : m *2;
                bps = realloc(bps, sizeof(int)*m);
            }
            if (bed0->end < bed->end) {
                bed0->end = bed->end;
                // push bed0->end
                bps[n++] = bed0->end;
            } else if (bed0->end > bed->end){
                //push bed->end
                bps[n++] = bed->end;
            }
            
            bed->seqname = -1;
        }
        if (n > 0) {
            if (n > 1) qsort(bps, n, sizeof(int), cmpint);
            int k;
            int last_start = bed0->start;
            for (k = 0; k < n; ++k) {
                bed_spec_push1(new, bed0->seqname, last_start-offset, bps[k], bed0->strand, bed0->name, NULL);
                last_start = bps[k];
            }
            bed_spec_push1(new, bed0->seqname, last_start-offset, bed0->end, bed0->strand, bed0->name, NULL);
        } else {
            bed_spec_push1(new, bed0->seqname, bed0->start-offset, bed0->end, bed0->strand, bed0->name, NULL);
        }
    }

    if (m>0) free(bps);
    
    bed_spec_destroy(B);
    bed_spec_dedup(new, 1);
    return new;
}
// level: 1 for gene, 2 for transcript, 3 for exon, 4 for CDS
// name_level : 0 for none, 1 for gene, 2 for transcript name, 3 for exon/CDS name
struct bed_spec *gtf2bed(struct gtf_spec *G, int level, int name_level, int offset)
{
    struct bed_spec *B = bed_spec_init();
    bed_spec_seqname_from_gtf(B, G);

    if (name_level > level) error("name_level should <= level.");

    if (offset != 0) offset = 1;
    
    int i;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        char *seqname = dict_name(G->name, i);
        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *g = ctg->gtf[j];

            char *name = NULL;
            
            if (name_level == 1) name = dict_name(G->gene_name, g->gene_name);
            
            if (level == 1) {
                bed_spec_push0(B, seqname, g->start-offset, g->end, g->strand, name, NULL);            
            } else {
                int k;
                for (k = 0; k < g->n_gtf; ++k) {
                    struct gtf *t = g->gtf[k];
                    if (name_level == 2) name = dict_name(G->transcript_id, t->transcript_id);
                    else if (name_level == 3) name = NULL;

                    if (level== 2) {
                        bed_spec_push0(B, seqname, t->start-offset, t->end, t->strand, name, NULL);   
                    } else {
                        int e;
                        for (e = 0; e < t->n_gtf; ++e) {
                            struct gtf *ex= t->gtf[e];
                            if (level == 3 && ex->type == feature_exon) {
                                bed_spec_push0(B, seqname, ex->start-offset, ex->end, ex->strand, name, NULL);
                            } else if (level == 4 && ex->type == feature_CDS) {
                                bed_spec_push0(B, seqname, ex->start-offset, ex->end, ex->strand, name, NULL);
                            }
                        }
                    }
                }
            }
        }
    }
    return B;
}
// Bed Annotation
// sort by type and distance
static int cmpfunc2(const void *_a, const void *_b)
{
    const struct anno0 *a = (const struct anno0*) _a;
    const struct anno0 *b = (const struct anno0*) _b;
    if (a->type != b->type) return (a->type > b->type) - (a->type < b->type);
    return 0;
    // return (a->dist < b->dist) - (a->dist > b->dist);
}

// exon
static int query_exon(int start, int end, struct gtf const *G, struct anno0 *a, int coding)
{
    int utr = 0;
    if (coding) utr = 1; // if CDS record exists, turn utr to 0 if region overlapped with CDS region
    int pass_cds = 0;

    int n = 0;
    struct gtf **gtf_pool = malloc(G->n_gtf *sizeof(struct gtf*));
    int i;
    for (i = 0; i < G->n_gtf; ++i) { // exon level
        struct gtf *g0 = G->gtf[i];
        // filter type
        if (g0->type != feature_CDS && g0->type != feature_exon) continue;
        // check overlapped
        
        // non-overlap
        if (end <= g0->start) {
            if (n == 0) { // put this record, as intron
                gtf_pool[n++] = G->gtf[i];
            }
            break;   
        }

        if (g0->type == feature_CDS) {
            if (start > g0->end) pass_cds = 1;
        }
            
        if (start > g0->end) continue; // check next
        
        // CDS record is only used to distiguish UTR and EXON
        if (g0->type == feature_CDS) {
            utr = 0; // it's a coding region
            continue;
        }
        
        // push to pool
        gtf_pool[n++] = g0;

        // out of range, no need to check next one
        if (end <= g0->end) break;
    }
    
    if (n == 0) {
        a->type = BAT_INTERGENIC;            
        // should not come here
        a->g = NULL;
    }
    else if (n == 1) {
        struct gtf *g0 = gtf_pool[0];
        if (g0->start <= start && g0->end >= end) {
            a->type = BAT_EXON;
        }
        else if (end <= g0->start) {
            a->type = BAT_INTRON;
        }
        else {
            a->type = BAT_EXONINTRON;
        }
        a->g = g0;
    }
    else {
        struct gtf *g0 = gtf_pool[0];
        a->type = BAT_MULTIEXONS;
        a->g = g0;
    }

    if (a->type == BAT_EXON && utr == 1) {
        // forward
        if (G->strand == 0) a->type = pass_cds ? BAT_UTR3 : BAT_UTR5;
        // backward
        else a->type = pass_cds ? BAT_UTR5 : BAT_UTR3;
    }

    free(gtf_pool);
    
    return a->type == BAT_INTERGENIC;
}

static int query_trans(int start, int end, struct gtf const *G, struct anno0 *a)
{
    // nonoverlap with this gene
    if (start >= G->end) return BAT_FLANK;
    if (end <= G->start) return BAT_FLANK;

    // partial overlap
    if (start < G->start) return BAT_FLANK;
    if (end > G->end) return BAT_FLANK;
    
    // debug_print("gene: %s, %d\t%d", GTF_genename(args.G,G->gene_name), G->start, G->end);
    struct anno0 *a0 = malloc(sizeof(struct anno0)* G->n_gtf);
    int i;
    int j = 0;
    for (i = 0; i < G->n_gtf; ++i) { // 
        struct gtf *g0 = G->gtf[i];
        if (g0->type != feature_transcript) continue;
        
        //  nonoverlap with this trans
        if (start >= g0->end) continue;
        if (end <= g0->start) continue;
        
        int ret = query_exon(start, end, g0, &a0[j], g0->coding);
        if (ret == 0) j++;
    }
    
    if (j ==0) { free(a0); return -1; }
    // debug_print("j : %d",j);
    
    if (j > 1) qsort(a0, j, sizeof(struct anno0), cmpfunc2);

    // if hit more than one exon/cds record, assign by the first record
    a->type = a0[0].type;
    a->g = a0[0].g;
    // debug_print("%s\t%s", bed_typename(a->type), GTF_transid(args.G, a->g->transcript_id));
    free(a0);
    //debug_print("type: %d", a->type);
    return a->type;
}

/* static int query_promoter(int start, int end, struct gtf *G, struct anno0 *a, int down, int up) */
/* { */
/*     // nonoverlap with this promoter */
/*     int start0 = G->start; */
/*     int end0 = G->start; */
/*     if (G->strand == 1) { */
/*         start0 = G->end; */
/*         end0 = G->end; */
/*         start0 = start0 - down; */
/*         end0 = start0 + up; */
/*     } else { */
/*         start0 = start0 - up; */
/*         end0 = start0 + down; */
/*     } */

/*     if (start0 < 0) start0 = 1; */
/*     if (end0 < 0) end0 = start0; */
    
/*     if (start >= end0) return -1; */
/*     if (end <= start0) return -1; */

/*     a->type = BAT_PROMOTER; */
/*     a->g = G; */
/*     return BAT_PROMOTER; */
/* } */

#define max(x, y) x > y ? x : y
#define min(x, y) x > y ? y : x

int region_overlap(int start, int end, int start0, int end0)
{
    int start1 = min(start, start0);
    int end1 = max(end, end0);

    return end - start + end0 - start0 - (end1 -start1);
}

static struct dict *wnames = NULL;

void anno_bed_cleanup()
{
    if (wnames !=NULL) dict_destroy(wnames);    
}
struct anno0 *anno_bed_core(const char *name, int start, int end, int strand, struct gtf_spec *G, int *n, int down, int up)
{
    if (end == -1) end = start+1;
    if (end < start ) {
        warnings("end < start");
        return NULL;
    }

    int flank = max(down, up);
    
    *n = 0;
    struct region_itr *itr = gtf_query(G, name, start - flank, end + flank);
    if (itr == NULL) {
        int id = dict_query(G->name, name);
        if (id == -1) {
            if (wnames == NULL) wnames = dict_init();
            int idx = dict_query(wnames, name);
            if (idx == -1) {
                warnings("Chromosome %s not found in GTF, use wrong database? ", name);
                dict_push(wnames, name);
            }
        }
        return NULL;
    }
    
    // annotate all possibility
    struct anno0 *a = malloc(sizeof(struct anno0)*itr->n);
    int k = 0;
    for (int j = 0; j < itr->n; ++j) {
        struct gtf *g0 = (struct gtf*)itr->rets[j];
        if (start <= g0->start && end >= g0->end) {
            a[k].type = BAT_WHOLEGENE;
            a[k].g = g0;
        } else {
            int ret;
            /* if (promoter) { */
            /*     ret = query_promoter(start, end, g0, &a[k], down, up); */
            /*     if (ret == -1) { */
            /*         ret = query_trans(start, end, g0, &a[k]); */
            /*     } */
            /* } else { */
            ret = query_trans(start, end, g0, &a[k]);
                //}

            if (ret == BAT_FLANK) {
                /* if (strand == -1) { */
                /*     a[k].type = BAT_INTERGENIC; */
                /*     continue; */
                /* }  */

                /* if (g0->strand == strand) { */
                /*     a[k].type = BAT_INTERGENIC; */
                /*     continue; */
                /* } */

                if (strand_is_minus(g0->strand)) {
                    if (start >= g0->end + up) {
                        a[k].type = BAT_INTERGENIC;
                        continue;
                    }

                    if (end <= g0->start - down) {
                        a[k].type = BAT_INTERGENIC;
                        continue;
                    }

                    
                    if (region_overlap(start, end, g0->end, g0->end + up) > 0) {
                        a[k].type = BAT_UPSTREAM;
                        a[k].g = g0;
                    } else if (region_overlap(start, end, g0->start - down, g0->start) > 0) {
                        a[k].type = BAT_DOWNSTREAM;
                        a[k].g = g0;
                    }
                } else {
                    if (start >= g0->end + down) {
                        a[k].type = BAT_INTERGENIC;
                        continue;
                    }
                    if (end <= g0->start - up) {
                        a[k].type = BAT_INTERGENIC;
                        continue;
                    }

                    if (region_overlap(start, end, g0->start - up, g0->start) > 0) {
                        a[k].type= BAT_UPSTREAM;
                        a[k].g = g0;
                    } else if (region_overlap(start, end, g0->end, g0->end + down) >0) {
                        a[k].type = BAT_DOWNSTREAM;
                        a[k].g = g0;
                    }
                }
            }
            
            if (ret == -1) continue; // nonoverlapped
        }

        // for stranded, reannotate antisense 
        if (strand != -1 && g0->strand != strand) {
            if (a[k].type == BAT_MULTIEXONS) a[k].type = BAT_ANTISENSECOMPLEX;
            else if (a[k].type == BAT_WHOLEGENE) a[k].type = BAT_ANTISENSECOMPLEX;
            else if (a[k].type == BAT_EXONINTRON) a[k].type = BAT_ANTISENSECOMPLEX;
            else if (a[k].type == BAT_UTR3) a[k].type = BAT_ANTISENSEUTR3;
            else if (a[k].type == BAT_UTR5) a[k].type = BAT_ANTISENSEUTR5;
            else if (a[k].type == BAT_EXON) a[k].type = BAT_ANTISENSEEXON;
            else if (a[k].type == BAT_INTRON) a[k].type = BAT_ANTISENSEINTRON;
            else if (a[k].type == BAT_UPSTREAM) a[k].type = BAT_ANTISENSEUP;
            else if (a[k].type == BAT_DOWNSTREAM) a[k].type = BAT_ANTISENSEDOWN;
        }

        k++;
    }

    region_itr_destroy(itr);

    if (k > 1) {
        qsort(a, k, sizeof(struct anno0), cmpfunc2);
        for (int j = 1; j < k; ++j) { //incase multiple antisense situations
            if (a[j].type > BAT_WHOLEGENE) {
                k = j;
                break;
            }
        }
    }

    if (k == 0) {
        free(a);
        return NULL;
    }
    
    *n = k;
    return a;
}
