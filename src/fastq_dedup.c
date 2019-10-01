#include "utils.h"
#include "fastq.h"
#include "htslib/thread_pool.h"

const char *dedup_version = "0.0.0.9999";

static int usage()
{
    fprintf(stderr, "Programs : Deduplicate reads with same barcodes.\n");
    fprintf(stderr, "Version : %s\n", dedup_version);
    fprintf(stderr, "Options :\n");
    fprintf(stderr, "  -tag    Tags of each block, usually be cell barcode and UMI. Such as CB,UY.");
    fprintf(stderr, "  -o      Output reads.\n");
    fprintf(stderr, "  -t      Threads.\n");
    return 1;
}

static void *run_it(void *_d)
{
    
}
                        
int fastq_dedup(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();
    
            hts_tpool *pool = hts_tpool_init(args.n_thread);
        hts_tpool_process *q = hts_tpool_process_init(pool, args.n_thread*2, 0);
        hts_tpool_result *r;        
        for (i = 0; i < idx->n; ++i) {
            char *key = idx->key[i];        
            struct bseq_pool *p = bseq_pool_init();        
            
            k = kh_get(idx,(kh_idx_t*)idx->dict, key);
            struct data_index *di = kh_val((kh_idx_t*)idx->dict, k);
            for (j = 0; j < di->n; ++j) {
                struct bseq *b = read_file_block(args.fp_in, &di->idx[j]);
                bseq_pool_push(b, p);
            }
            int block;
            do {
                block = hts_tpool_dispatch2(pool, q, run_it1, p, 1);
                if ((r = hts_tpool_next_result(q))) {                
                    struct bseq_pool *d = (struct bseq_pool*)hts_tpool_result_data(r);
                    write_out(d);
                    hts_tpool_delete_result(r, 0);
                }
            }
            while (block == -1);
        }
        hts_tpool_process_flush(q);
        
        while ((r = hts_tpool_next_result(q))) {
            struct bseq_pool *d = (struct bseq_pool*)hts_tpool_result_data(r);
            write_out(d);
            hts_tpool_delete_result(r, 0);
        }
        hts_tpool_process_destroy(q);
        hts_tpool_destroy(pool);

}
