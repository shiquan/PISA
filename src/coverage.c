#include "utils.h"
#include "coverage.h"
#include "bed.h"

struct depth *depth_init()
{
    struct depth *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    d->id = -1;
    return d;
}

void update_depth(struct depth *d)
{
    if (d->bc1 && d->dep1 == 0) {
        d->dep1 = dict_size(d->bc1);
        dict_destroy(d->bc1);
    }
    if (d->bc2 && d->dep2 == 0) {
        d->dep2 = dict_size(d->bc2);
        dict_destroy(d->bc2);
    }
}

int depth_increase(struct depth *d, int id, char *umi, int strand0)
{
    if (d->id == -1) d->id = id;
    if (d->id == id) {
        if (strand0 == BED_STRAND_UNK || strand0 == BED_STRAND_FWD) {
            if (umi) {
                if (d->bc1 == NULL) d->bc1 = dict_init();
                dict_push(d->bc1, umi);
            }
            else d->dep1++;
        } else {
            if (umi) {
                if (d->bc2 == NULL) d->bc2 = dict_init();
                dict_push(d->bc2, umi);
            }
            else d->dep2++;
        }
        
    } else if (d->id > id) {
        if (d->right) return depth_increase(d->right, id, umi, strand0);
        int pos = d->pos;
        struct depth *new = depth_init();
        new->pos = pos;
        new->id = id;
        if (strand0 == BED_STRAND_UNK || strand0 == BED_STRAND_FWD) {
            if (umi) {
                if (new->bc1 == NULL) new->bc1 = dict_init();
                dict_push(new->bc1, umi);
            }
            else new->dep1++;
        } else {
            if (umi) {
                if (new->bc2 == NULL) new->bc2 = dict_init();
                dict_push(new->bc2, umi);
            }
            else new->dep2++;
        }
        d->right = new;
    } else { // cur->id < id
        if (d->left) return depth_increase(d->left, id, umi, strand0);
        int pos = d->pos;
        struct depth *new = depth_init();
        new->pos = pos;
        new->id = id;
        if (strand0 == BED_STRAND_UNK || strand0 == BED_STRAND_FWD) {
            if (umi) {
                if (new->bc1 == NULL) new->bc1 = dict_init();
                dict_push(new->bc1, umi);
            }
            else new->dep1++;
        } else {
            if (umi) {
                if (new->bc2 == NULL) new->bc2 = dict_init();
                dict_push(new->bc2, umi);
            }
            else new->dep2++;
        }
        d->left = new;
    }
    return 0;
}

struct depth* bam2depth(const hts_idx_t *idx, const int tid, const int start, const int end,
                        const int strand,
                        htsFile *fp,
                        const int mapq_thres,
                        const int ignore_strand,
                        struct dict *bc,
                        const char *bc_tag,
                        const char *umi_tag,
                        const int split_by_tag,
                        const int alias_tag,
                        const int *alias_idx
    )
{
    hts_itr_t *itr = sam_itr_queryi(idx, tid, start, end);
    if (itr == NULL) return NULL;
    bam1_t *b = bam_init1();

    int result;

    struct depth *head = NULL;
    struct depth *tail = NULL;

    struct depth *last_start = NULL; // point to start pos of last read
    //struct depth *cur0 = NULL;
    int strand0;
    
    // int last = start+1;
        
    while ((result = sam_itr_multi_next(fp, itr, b)) >= 0) {
        bam1_core_t *c = &b->core;
        if (c->qual < mapq_thres) continue;
        if (c->flag & BAM_FSECONDARY) continue;
        if (c->flag & BAM_FUNMAP) continue;
        if (c->flag & BAM_FQCFAIL) continue;
        if (c->flag & BAM_FDUP) continue;
        
        if (ignore_strand) strand0 = BED_STRAND_UNK;
        else {
            if (c->flag & BAM_FREVERSE) strand0 = BED_STRAND_REV;
            else strand0 = BED_STRAND_FWD;
        }
        
        if (strand == BED_STRAND_FWD && strand0 == BED_STRAND_REV) continue;
        if (strand == BED_STRAND_REV && strand0 == BED_STRAND_FWD) continue;
        
        uint8_t *data = NULL;
        char *umi = NULL;
        int id = -1;
        
        if (bc_tag && bc) {
            data = bam_aux_get(b, bc_tag);
            if (data == NULL) continue;
            id = dict_query(bc, (char*)(data+1));
            if (id == -1) continue;
            if (alias_tag) {
                id = alias_idx[id];
            }
        }
        if (split_by_tag) {
            if (id == -1) {
                data = bam_aux_get(b, bc_tag);
                if (data == NULL) continue;
                id = dict_query(bc, (char*)(data+1));
                if (id == -1) id = dict_push(bc, (char*)(data+1));
            }
        } else {
            id = -1;
        }
        
        if (umi_tag) {
            data = bam_aux_get(b, umi_tag);
            if (data == NULL) continue;
            umi = (char*)(data+1);
        }
        
        int pos = c->pos+1;
        struct depth *cur = NULL;
        int i;
        for (i = 0; i < b->core.n_cigar; ++i) {
            int cig = bam_cigar_op(bam_get_cigar(b)[i]);
            int ncig = bam_cigar_oplen(bam_get_cigar(b)[i]);
            if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
                int j;
                for (j = 0; j < ncig; ++j) {
                    if (pos <= start) { pos++; continue;} // overlapped but not fully enclosed reads, at this moment last_start should be NULL
                    if (pos > end) { i = b->core.n_cigar; break; } // skip this record
                    
                    if (head == NULL) {
                        head = depth_init();
                        head->pos = pos;
                        tail = head;
                        // head is not always point to last_start, consider spliced reads may start from far away positions
                        //last_start = head;
                    }

                    // this is for fully enclosed or left covered reads 
                    if (cur && cur->pos == c->pos+1) last_start = cur; // start pos of current record, start to compare for next record
                    
                    //         header     tail
                    //         |          |
                    // id:     123412341234  for cell id ...
                    // pos:    111122223333  for each position ...
                    //
                    if (cur == NULL) { // first pos, track from start pos of last read
                        // debug_print("new : %s", bam_get_qname(b));
                        cur = last_start == NULL ? head : last_start; // if right covered, start from head for safty
                        // debug_print("last: %d, pos: %d", cur->pos, pos);
                        // assert(cur->pos <= pos);
                        for (; cur != NULL; cur = cur->next) {
                            if (cur->pos == pos) {                                
                                depth_increase(cur, id, umi, strand0);
                                break;
                            }
                            else if (cur->pos > pos) { // spliced pos
                                struct depth *new = depth_init();
                                new->pos = pos;
                                new->next = cur;
                                new->before = cur->before;
                                if (cur->before) cur->before->next = new;
                                else {
                                    assert(cur == head);
                                    head->before = new;
                                    head = new;
                                }
                                cur->before = new;
                                cur = new;
                                depth_increase(cur, id, umi, strand0);
                                break;
                            }
                        }
                        
                        if (cur == NULL) { // update tail
                            struct depth *new = depth_init();
                            new->pos = pos;                         
                            new->before = tail;
                            tail->next = new;
                            tail = new;
                            cur = new;
                            depth_increase(cur, id, umi, strand0);
                        }
                    }
                    else { // not the first one, count start from cur
                        for (;;) {
                            if (cur == NULL) break;
                            if (cur->pos == pos) {
                                depth_increase(cur, id, umi, strand0);
                                break;
                            }
                            else if (cur->pos > pos) {
                                struct depth *new = depth_init();
                                new->pos = pos;
                                new->next = cur;
                                new->before = cur->before;
                                if (cur->before) cur->before->next = new;
                                else {
                                    //debug_print("cur: %d, head: %d", cur->pos, head->pos);
                                    assert(cur == head);
                                    head->before = new;
                                    head = new;
                                }
                                cur->before = new;
                                cur = new;
                                depth_increase(cur, id, umi, strand0);
                                break;
                            }
                            cur = cur->next; // next record
                        }
                            
                        if (cur == NULL) {
                            assert(tail);
                            struct depth *new = depth_init();
                            new->pos = pos;
                            new->before = tail;
                            tail->next = new;
                            tail = new;
                            cur = new;
                            depth_increase(cur, id, umi, strand0);
                        }
                    }
                    pos++;
                }
            }
            else if (cig == BAM_CDEL) {
                pos += ncig;
            }
            else if (cig == BAM_CREF_SKIP) {
                pos += ncig;
            }
        }

        // assert(tail);
        // debug_print("%s\t%d", bam_get_qname(b), last_start ? last_start->pos : -1);
    }
    
    if (result < -1)
        warnings("Failed to retrieve region due to truncated file or corrupt bam index.");
    
    hts_itr_destroy(itr);

    struct depth *cur = head;
    while (cur) {
        update_depth(cur);
        cur = cur->next;
    }
    
    bam_destroy1(b);
    return head;
}

