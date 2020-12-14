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

int reads_match_var(struct bed *bed, bam1_t *b)
{
    int st = bam_seqpos(b, bed->start);
    if (st == -1) return -1; // out of range
    if (st == -2) return -1; // intron region

    struct var *v = (struct var*)bed->data;
    assert(v);

    if (st == -3) { // deletion
        if (v->alt->l < v->ref->l) return 0; // matched
        return 1; // unmatched
    } 
    uint8_t *seq = bam_get_seq(b);
    bam1_core_t *c = &b->core;
    int i, j;
    int mis = 0;
    for (i = st, j = 0; i < c->l_qseq && j < v->ref->l; ++i,++j) {
        if ("=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)] != v->ref->s[j]) {
            mis = 1;
            break;
        }
    }

    if (mis) {
        for (i = st, j = 0; i < c->l_qseq && j < v->alt->l; ++i,++j) {
            if ("=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)] != v->alt->s[j])
                return 1;
        }
        return 0;
    }
            
    return 1;       
}
int bam_vcf_anno(bam1_t *b, bam_hdr_t *h, struct bed_spec const *B, const char *vtag)
{ 
    bam1_core_t *c;
    c = &b->core;
    
    // cleanup exist tag
    uint8_t *data;
    if ((data = bam_aux_get(b, vtag)) != NULL) bam_aux_del(b, data);
    
    char *name = h->target_name[c->tid];
    int endpos = bam_endpos(b);

    struct region_itr *itr = bed_query(B, name, c->pos, endpos);
    if (itr == NULL) return 0; // query failed
    if (itr->n == 0) return 0; // no hit

    struct dict *val = dict_init();
    int i;
    kstring_t temp = {0,0,0};
    for (i = 0; i < itr->n; ++i) {
        struct bed *bed = (struct bed*)itr->rets[i];
        if (bed->start > endpos || bed->end <= c->pos) continue; // not covered
        if (reads_match_var(bed, b) != 0) continue; // no variant contained
        
        temp.l = 0;

        if (bed->name == -1) ksprintf(&temp, "%s,%d,%s", dict_name(B->seqname, bed->seqname), bed->start+1, ((struct var*)bed->data)->alt->s);
        else kputs(dict_name(B->name,bed->name), &temp);
        
        if (temp.l)
            dict_push(val, temp.s);
        
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
