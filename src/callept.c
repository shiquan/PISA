#include "utils.h"
#include "htslib/sam.h"
#include "coverage.h"
#include "bed.h"
#include "number.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bc_list;
    const char *tag;
    const char *umi_tag;

    struct dict *barcodes;

    htsFile *fp;
    hts_idx_t *idx;
    FILE *out;

    // struct bed_spec *B;

    bam_hdr_t *hdr;

    int mapq_thres;
    int n_thread;

    int ignore_strand;

    int min_length;
    int max_gap;
    int cutoff;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .bc_list = NULL,
    .tag = NULL,
    .umi_tag = NULL,
    .barcodes = NULL,
    .fp = NULL,
    .idx = NULL,
    .out = NULL,
    // .B = NULL,
    .hdr = NULL,
    .mapq_thres = 20,
    .n_thread = 4,
    .ignore_strand = 0,
    .min_length = 50,
    .max_gap = 50,
    .cutoff = 1
};

static int callept_usage()
{
    fprintf(stderr, "# Call Expressed Peak Tags (EPTs) for indexed BAMs.\n");
    fprintf(stderr, "# Report peaks which depths are continuously higher than a given cutoff.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m cellept -o epts.bed sorted.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m cellept -tag CB -list cells.txt -umi UB -o epts.bed sorted.bam\n");
    fprintf(stderr, "\nOptions : \n");
    fprintf(stderr, " -tag      [TAG]      Tag used for grouping reads.\n");
    fprintf(stderr, " -list     [FILE]     Candidate list for -tag.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. If set, only count unique UMIs for each location.\n");
    fprintf(stderr, " -is                  Ignore strand.\n");
    fprintf(stderr, " -gap      [INT]      Maximum gap to merge nearby peaks. [50]\n");
    fprintf(stderr, " -min-length  [INT]   Minimum peak length. [200]\n");
    fprintf(stderr, " -cutoff   [INT]      Cutoff of depth. [1]\n");
    fprintf(stderr, " -o        [FILE]     Output EPTs in bed format. [stdout].\n");
    fprintf(stderr, " -q        [INT]      Minimal map quality to filter. [20]\n");
    fprintf(stderr, " -t        [INT]      Threads. [4]\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * Requires sorted and indexed BAM as input.\n");
    fprintf(stderr, " * Compares with `MACS2` and other peak callers, PISA callept considers UMIs and strand of reads.\n");
    return 1;
}

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    const char *threads = NULL;
    const char *cutoff = NULL;
    const char *gap = NULL;
    const char *minlength = NULL;
    
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-list") == 0) var = &args.bc_list;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-t") == 0) var = &threads;        
        else if (strcmp(a, "-is") == 0) {
            args.ignore_strand = 1;
            continue;
        }
        else if (strcmp(a, "-gap") == 0) var = &gap;
        else if (strcmp(a, "-min-length") == 0) var = &minlength;
        else if (strcmp(a, "-cutoff") == 0) var = &cutoff;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument, %s", a);
    }

    if (args.input_fname == 0) error("No input bam.");

    if (threads) args.n_thread = str2int((char*)threads);
    if (gap) args.max_gap = str2int((char*)gap);
    if (minlength) args.min_length = str2int((char*)minlength);
    if (cutoff) args.cutoff = str2int((char*)cutoff);
    
    if (mapq) {
        args.mapq_thres = str2int(mapq); 
    }

    if (args.tag) args.barcodes = dict_init();
    
    if (args.bc_list) {
        int val = 0;
        dict_read2(args.barcodes, args.bc_list, &val);
        if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");
    }
    
    args.fp = hts_open(args.input_fname, "r");
    if (args.fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    
    hts_set_threads(args.fp, args.n_thread);

    args.idx = sam_index_load(args.fp, args.input_fname);
    if (args.idx == NULL) error("Failed to load bam index of %s", args.input_fname);

    args.hdr = sam_hdr_read(args.fp);
    if (args.output_fname) {
        args.out = fopen(args.output_fname, "w");
        if (args.out == NULL) error("%s : %s.", args.output_fname, strerror(errno));        
    }
    else args.out = stdout;

    return 0;
}

static void memory_release()
{
    if (args.barcodes) dict_destroy(args.barcodes);
    
    bam_hdr_destroy(args.hdr);
    hts_idx_destroy(args.idx);
    hts_close(args.fp);
    // if (args.B) bed_spec_destroy(args.B);
    
    if (args.out != stdout) fclose(args.out);    
}

const static int ws = 100000;

int callept(struct bed_spec *B, hts_idx_t *idx, int tid, int start, int end)
{
    struct depth *d = bam2depth(idx, tid, start, end, BED_STRAND_UNK, args.fp, args.mapq_thres,
                                args.ignore_strand, args.barcodes, args.tag, args.umi_tag,
                                0, 0, NULL);

    if (d == NULL) return 1;
    
    struct bed b1 = {tid, -1, 0, 0, 0, NULL};
    struct bed b2 = {tid, -1, 0, 0, 1, NULL};

    for (;;) {
        if (d == NULL) break;
        
        if (d->dep1 > args.cutoff || d->dep2 > args.cutoff) {
            //debug_print("%d\n", d->pos);
            assert(d->pos > b1.end);
            assert(d->pos > b2.end);
            
            // update new region
            if (d->pos > b1.end + 1 && b1.end > 0) {
                bed_spec_push(B, &b1);
                b1.start = 0;
                b1.end = 0;
                //b1.start = d->pos-1; // 0 based
                //b1.end = d->pos;
            }

            if (d->pos > b2.end + 1 && b2.end > 0) {
                bed_spec_push(B, &b2);
                b2.start = 0;
                b2.end = 0;
                //b2.start = d->pos-1; // 0 based
                //b2.end = d->pos;
            }
            
            // span region
            if (d->dep1 < args.cutoff) {
                if (b1.end > 0) {
                    bed_spec_push(B, &b1);
                    b1.start = 0;
                    b1.end = 0;
                }
            } else {
                if (b1.start == 0 && b1.end == 0) {
                    b1.start = d->pos -1;
                    b1.end = d->pos;
                } else {
                    assert (d->pos == b1.end + 1);
                    b1.end ++; // span one
                }
            }

            if (d->dep2 < args.cutoff) {
                if (b2.end > 0) {
                    bed_spec_push(B, &b2);
                    b2.start = 0;
                    b2.end = 0;
                }
            } else {
                if (b2.start == 0 && b2.end == 0) {
                    b2.start = d->pos -1;
                    b2.end = d->pos;
                } else {
                    assert (d->pos == b2.end + 1);
                    b2.end ++;
                }
            }
            // move to next pos
            // continue;
        }

        struct depth *tmp = d;
        d = d->next;
        free(tmp);
    }

    if (b1.end > 0) bed_spec_push(B, &b1);
    if (b2.end > 0) bed_spec_push(B, &b2);
    
    return 0;
}
int callept_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return callept_usage();
    struct bed_spec *B = bed_spec_init();
    bed_spec_seqname_from_bam(B, args.hdr);
    
    int i;
    for (i = 0; i < args.hdr->n_targets; ++i) {
        LOG_print("Process %s ..", args.hdr->target_name[i]);
        int len = args.hdr->target_len[i];
        int last = 0;
        int end = last + ws;
        if (end > len) end = len;
        
        for (;;) {
            callept(B, args.idx, i, last, end);

            if (end == len) break;
            last = end;
            end = end + ws;
            if (end > len) end = len;
        }
    }

    bed_spec_merge2(B, 1, args.max_gap, args.min_length);
    bed_spec_write0(B, args.out);

    bed_spec_destroy(B);
    memory_release();
    
    return 0;
}
