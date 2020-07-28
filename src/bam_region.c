#include "bam_region.h"
#include "bed.h"
struct bam_region_itr *bri_init(htsFile *fp, const char *fn, struct bed_spec *bed)
{
    if (bed == NULL) return NULL;
    struct bam_region_itr *bri = malloc(sizeof(*bri));
    memset(bri, 0, sizeof(*bri));
    bri->bed = bed;
    bri->idx = sam_index_load(fp, fn);
    if (bri->idx == NULL) error("Failed to load index of %s.", fn);
    return bri;
}

int bri_read(htsFile *fp, bam_hdr_t *hdr, bam1_t *b, struct bam_region_itr *bri)
{
    for (;;) {
        if (bri->i_bed == bri->bed->n) return -1; // end
        struct bed *bed = &bri->bed->bed[bri->i_bed++];
        int tid = sam_hdr_name2tid(hdr, dict_name(bri->bed->seqname, bed->seqname));
        bri->itr = sam_itr_queryi(bri->idx, tid, bed->start, bed->end);
        if (bri->itr) {
            int ret;
            if ((ret = sam_itr_next(fp, bri->itr, b)) >= 0) return ret;
            hts_itr_destroy(bri->itr);
            bri->i_bed++;
        }
    }
}

void bri_destroy(struct bam_region_itr *bri)
{
    bed_spec_destroy(bri->bed);
    hts_idx_destroy(bri->idx);
    free(bri);
}
