// pick alignment reads by cell barcodes
#include "utils.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include "barcode_list.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/khash_str2int.h"
#include "htslib/kseq.h"
#include "bam_pool.h"

static int usage()
{
    fprintf(stderr, "PickBam -list barcode.txt -tag CB -o out.bam in.bam\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *barcode_fname;
    const char *tag;
    int file_th;
    int n_thread;
    struct barcode_list *barcode;

    htsFile *in;
    htsFile *out;
    bam_hdr_t *hdr;

    int chunk_size;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .barcode_fname = NULL,
    .tag = NULL,
    .file_th = 4,
    .n_thread = 1,
    .barcode = NULL,

    .in = NULL,
    .out = NULL,
    .hdr = NULL,
    .chunk_size = 1000000, // 1M
};
static int parse_args(int argc, char **argv)
{
    const char *file_th = NULL;
    const char *thread = NULL;
    int i;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-list") == 0) var = &args.barcode_fname;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-t") == 0) var = &thread;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument: %s", a);
    }

    // CHECK_EMPTY(args.barcode_fname, "-list Cell barcode list must be set.");
    CHECK_EMPTY(args.input_fname, "Input BAM file must be set.");
    CHECK_EMPTY(args.output_fname, "Output BAM file must be set.");

    if (file_th) args.file_th = str2int((char*)file_th);
    if (thread) args.n_thread = str2int((char*)thread);
    
    if (args.barcode_fname) {
        args.barcode = barcode_init();
        if (barcode_read(args.barcode, args.barcode_fname)) error("Empty barcode");
    }

    args.in = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.in, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    hts_set_threads(args.in, args.file_th);
    //hts_set_threads(args.out, args.file_th);

    args.hdr = sam_hdr_read(args.in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");

    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
    return 0;
}
static void memory_release()
{
    if (args.barcode)
        barcode_destory(args.barcode);
    bam_hdr_destroy(args.hdr);
    sam_close(args.in);
    sam_close(args.out);
}
static void write_out(struct bam_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = &p->bam[i];
        if (b->core.flag & BAM_FQCFAIL) continue;
        if (sam_write1(args.out, args.hdr, b) == -1) error("Failed to write SAM.");
    }
    bam_pool_destory(p);
}
static void *run_it(void *data)
{
    struct bam_pool *p = (struct bam_pool*)data;    
    int i;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = &p->bam[i];
        uint8_t *tag = bam_aux_get(b, args.tag);
        if (!tag) {
            b->core.flag = BAM_FQCFAIL;
            continue;
        }

        if (args.barcode) {
            int id;
            id = barcode_select(args.barcode, (char*)(tag+1));
            if (id == -1)
                b->core.flag = BAM_FQCFAIL; // destroy this bam
        }
    }
    return p;
}
    
int bam_pick(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return usage();
        
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    for (;;) {
        struct bam_pool *b = bam_pool_create();
        bam_read_pool(b, args.in, args.hdr, args.chunk_size);
            
        if (b == NULL) break;
        if (b->n == 0) { free(b->bam); free(b); break; }
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
                write_out(d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    memory_release();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());

    return 0;
}
