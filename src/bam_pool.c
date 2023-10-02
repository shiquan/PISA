#include "utils.h"
#include "bam_pool.h"

struct bam_pool *bam_pool_create()
{
    struct bam_pool *p = malloc(sizeof(*p));
    p->n = p->m =0;
    p->bam = NULL;
    return p;
}
struct bam_pool *bam_pool_init(int size)
{
    struct bam_pool *p = malloc(sizeof(*p));
    p->m = size;
    p->n = 0;
    
    if (p->m > 0) {
        p->bam = malloc(p->m*sizeof(bam1_t));
        memset(p->bam, 0, p->m*sizeof(bam1_t));
    } else {
        p->bam = NULL;    
    }
    
    return p;
}

void bam_read_pool(struct bam_pool *p, htsFile *fp, bam_hdr_t *h, int chunk_size)
{
    p->n = 0;
    int ret = -1;
    do {
        if (p->n >= chunk_size) break;
        if (p->n == p->m) {
            p->m = chunk_size;
            p->bam = realloc(p->bam, p->m*sizeof(bam1_t));
            for (int i = p->n; i <p->m; ++i) memset(&p->bam[i], 0, sizeof(bam1_t));
        }
        
        ret = sam_read1(fp, h, &p->bam[p->n]);
        if (ret < 0) break;
        p->n++;
    } while(1);
    //p->hdr = h;
    if (ret < -1) warnings("Truncated file?");    
}
void bam_pool_destory(struct bam_pool *p)
{
    int i;
    for (i = 0; i <p->n; ++i) 
        free(p->bam[i].data);
    free(p->bam);
    free(p);
}

