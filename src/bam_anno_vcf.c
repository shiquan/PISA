#include "utils.h"
#include "bed.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"

// return -1 on out of range
// -2 on intron region
// -3 on deletion
static int bam_seqpos(bam1_t *a, int pos)
{
    int st = 0;
    int i;
    const bam1_core_t *c = &a->core;
    int qpos = c->pos;
    uint32_t *cigar = bam_get_cigar(a);
    for (i = 0; i < c->n_cigar; ++i) {
        int type = cigar[i]&0xf;
        int len = cigar[i]>>4;
        if (type == BAM_CMATCH || type == BAM_CEQUAL || type == BAM_CDIFF) {
            int k;
            for (k = 0; k < len; ++k) {
                if (qpos == pos) return st;
                qpos++;
                st++;
            }
        }
        else if (type == BAM_CINS) {
            st+=len;
        }
        else if (type == BAM_CDEL) {
            qpos+=len;
            if (qpos > pos) return -3;
        }
        else if (type == BAM_CREF_SKIP) {
            qpos+=len;
            if (qpos > pos) return -2;
        }
        else {
            st+=len;
        }
        if (qpos > pos) break;
    }

    return -1;
}

/* // 0 on ref, 1 on mismatch, -1 on other situation */
/* int reads_match_var(struct bed *bed, bam1_t *b) */
/* { */
/*     int st = bam_seqpos(b, bed->start); */
/*     if (st == -1) return -1; // out of range */
/*     if (st == -2) return -1; // intron region */

/*     struct var *v = (struct var*)bed->data; */
/*     assert(v); */

/*     if (st == -3) { // deletion */
/*         if (v->alt->l < v->ref->l) return 0; // matched */
/*         return 1; // unmatched */
/*     }  */
/*     uint8_t *seq = bam_get_seq(b); */
/*     bam1_core_t *c = &b->core; */
/*     int i, j; */
/*     int mis = 0; */
/*     for (i = st, j = 0; i < c->l_qseq && j < v->ref->l; ++i,++j) { */
/*         if ("=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)] != v->ref->s[j]) { */
/*             mis = 1; */
/*             break; */
/*         } */
/*     } */

/*     if (mis == 0) return 0; */
/*     if (mis) { */
/*         for (i = st, j = 0; i < c->l_qseq && j < v->alt->l; ++i,++j) { */
/*             if ("=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)] != v->alt->s[j]) */
/*                 return -1; */
/*         } */
/*         return 1; */
/*     } */
            
/*     return 1;     */
/* } */

//  -1 on failed, 0 on ref, other number on allele index
int reads_match_var(struct bed *bed, bam1_t *b)
{
    int st = bam_seqpos(b, bed->start);
    if (st == -1) return -1; // out of range
    if (st == -2) return -1; // intron region

    //struct var *v = (struct var*)bed->data;
    bcf1_t *v = (bcf1_t*)bed->data;
    //bcf_hdr_t *hdr = (bcf_hdr_t*)bed->
    assert(v);
    bcf_unpack(v, BCF_UN_STR);

    if (st == -3) {
        int a;        
        for (a = 0; a < v->n_allele; ++a) {
            char *allele = v->d.allele[a];
            if (allele[0] == '*' && allele[1] == 0) return a;
        }
        return -1;
    }
        
    uint8_t *seq = bam_get_seq(b);
    bam1_core_t *c = &b->core;

    int a;        
    for (a = 0; a < v->n_allele; ++a) {
        int i, j;
        int mis = 0;
        char *allele = v->d.allele[a];
        for (i = st, j = 0; i < c->l_qseq && allele[j]; ++i,++j) {
            if ("=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)] != allele[j]) {
                mis = 1;
                break;
            }
        }
        if (mis == 0) return a;
    }

    return -1;
}

char *vcf_tag_name(int allele, bcf1_t *v, bcf_hdr_t *hdr, const struct bed_spec *B, struct bed *bed, int phased)
{
    kstring_t str={0,0,0};

    while (phased && v->n_sample > 0) {
        bcf_unpack(v, BCF_UN_FMT);
        
        bcf_fmt_t *fmt_ptr = bcf_get_fmt(hdr, v, "PGT");
        if (!fmt_ptr) break;
        
        int allele1=0, allele2=0;
        // int is_phased=0;
/*         int i; */

/* #define BRANCH_INT(type_t, vector_end) {                \ */
/*             type_t *p = (type_t*)(fmt_ptr->p);          \ */
/*             for (i=0; i<fmt_ptr->n; i++) {              \ */
/*                 if (p[i] == vector_end) break;          \ */
/*                 if (bcf_gt_is_missing(p[i])) break;     \ */
/*                 is_phased = p[i] &1;                    \ */
/*                 int tmp = p[i]>>1;                      \ */
/*                 if (tmp>1) {                            \ */
/*                     if (allele1==0) allele1=tmp;        \ */
/*                     else if (tmp != allele1) {          \ */
/*                         if (tmp < allele1) {            \ */
/*                             allele2 = allele1;          \ */
/*                             allele1 = tmp;              \ */
/*                         } else {                        \ */
/*                             allele2 = tmp;              \ */
/*                         }                               \ */
/*                     }                                   \ */
/*                 }                                       \ */
/*             }                                           \ */
/*         } */
        
/*         switch (fmt_ptr->type) { */
/*         case BCF_BT_INT8:  BRANCH_INT(int8_t, bcf_int8_vector_end); break; */
/*         case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break; */
/*         case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break; */
/*         default: error("Unexpected type %d", fmt_ptr->type); break; */
/*         } */
/*         #undef BRANCH_INT */
        char *src = (char*)(fmt_ptr->p);
        int l;
        for (l = 0; src[l] && l < fmt_ptr->size; ++l);
        if (l == fmt_ptr->size) break; // unknown phase block name
        assert(l>=3);

        kstring_t a1 = {0,0,0};
        kstring_t a2 = {0,0,0};

        int l0;
        for (l0 = 0; src[l0] != '|'; ++l0) kputc(src[l0], &a1);
        l0++;
        for (; src[l0] && l0 < l; ++l) kputc(src[l0], &a2);

        allele1 = str2int(a1.s);
        allele2 = str2int(a2.s);
        free(a1.s); free(a2.s);
        
        if (allele != allele1 && allele != allele2) break; // no phase allele exists
        fmt_ptr = bcf_get_fmt(hdr, v, "PID");
        if (!fmt_ptr) error("No PID tag found at phased position?");
        if (fmt_ptr->type != BCF_BT_CHAR) error("Inconsistant PID type.");
        src = (char*)(fmt_ptr->p);
        
        for (l = 0; src[l] && l < fmt_ptr->size; ++l);
        if (l == fmt_ptr->size) break; // unknown phase block name
        
        ksprintf(&str, "%s:%s-PB%d", dict_name(B->seqname, bed->seqname), src, allele);
        return str.s;   
    }
    
    ksprintf(&str, "%s:%d%s", dict_name(B->seqname, bed->seqname), bed->start+1, v->d.allele[0]);
    if (allele==0) {
        kputs("=", &str);        
    } else {
        kputc('>', &str);
        kputs(v->d.allele[allele],&str);
    }

    return str.s;
}
int bam_vcf_anno(bam1_t *b, bam_hdr_t *h, struct bed_spec const *B, const char *vtag, int ref_alt, int vcf_ss, int phased)
{ 
    bam1_core_t *c;
    c = &b->core;
    
    // cleanup exist tag
    uint8_t *data;
    if ((data = bam_aux_get(b, vtag)) != NULL) bam_aux_del(b, data);
    
    char *name = h->target_name[c->tid];
    int endpos = bam_endpos(b);
    
    // bed_query::start is 1 based
    struct region_itr *itr = bed_query(B, name, c->pos+1, endpos, BED_STRAND_IGN);
    if (itr == NULL) return 0; // query failed
    if (itr->n == 0) return 0; // no hit

    int strand = c->flag & BAM_FREVERSE;
    
    struct dict *val = dict_init();
    int i;
    kstring_t temp = {0,0,0};
    for (i = 0; i < itr->n; ++i) {
        struct bed *bed = (struct bed*)itr->rets[i];
        if (bed->start > endpos || bed->end <= c->pos) continue; // not covered
        int ret = reads_match_var(bed, b);
        temp.l = 0;

        if (ret == -1) continue;
        if (ref_alt == 1 && ret == 0) continue;
        
        char *name = vcf_tag_name(ret, (bcf1_t*)bed->data, (bcf_hdr_t*)B->ext, B, bed, phased);
        kputs(name, &temp);
        free(name);
        
        if (temp.l && vcf_ss) kputs(strand == 0 ? "/+" : "/-", &temp);
        if (temp.l) dict_push(val, temp.s);
    }
    
    region_itr_destroy(itr);

    if (temp.m) free(temp.s);
    
    if (dict_size(val)) {
        kstring_t str = {0,0,0};
        int i;
        for (i = 0; i < dict_size(val); ++i) {
            if (i) kputc(';', &str);
            kputs(dict_name(val, i), &str);
        }

        bam_aux_append(b, vtag, 'Z', str.l+1, (uint8_t*)str.s);
        free(str.s);
        dict_destroy(val);
        return 1;
    }
    dict_destroy(val);
    return 0;
}
