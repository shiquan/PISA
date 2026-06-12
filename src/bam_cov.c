// Calculate genomic coverage per cell/group tag
#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "dict.h"
#include "number.h"
#include "bed.h"
#include "coverage.h"

const static int ws1m = 1000000;

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bc_list;
    const char *tag;
    const char *umi_tag;

    struct dict *barcodes;
    int fix_barcodes;

    htsFile *fp;
    hts_idx_t *idx;
    bam_hdr_t *hdr;
    FILE *out;

    int mapq_thres;
    int min_depth;
    int n_thread;

    int is_csv;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .bc_list      = NULL,
    .tag          = NULL,
    .umi_tag      = NULL,
    .barcodes     = NULL,
    .fix_barcodes = 0,
    .fp           = NULL,
    .idx          = NULL,
    .hdr          = NULL,
    .out          = NULL,
    .mapq_thres   = 20,
    .min_depth    = 1,
    .n_thread     = 4,
    .is_csv       = 0,
};

extern int cov_usage();

static void write_output(uint64_t *covered, int n_cells, int sep_char)
{
    fprintf(args.out, "TAG%cCOVERED\n", sep_char);
    int i;
    for (i = 0; i < n_cells; ++i) {
        char *name = dict_name(args.barcodes, i);
        fprintf(args.out, "%s%c%" PRIu64 "\n", name, sep_char, covered[i]);
    }
}

static int parse_args(int argc, char **argv)
{
    int i;
    const char *mapq = NULL;
    const char *mindepth = NULL;
    const char *threads = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-list") == 0) var = &args.bc_list;
        else if (strcmp(a, "-umi") == 0) var = &args.umi_tag;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-q") == 0) var = &mapq;
        else if (strcmp(a, "-min-depth") == 0) var = &mindepth;
        else if (strcmp(a, "-@") == 0) var = &threads;

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
    if (args.tag == 0) error("No tag specified.");

    if (threads) args.n_thread = str2int((char*)threads);
    if (mapq) args.mapq_thres = str2int((char*)mapq);
    if (mindepth) args.min_depth = str2int((char*)mindepth);

    args.barcodes = dict_init();

    if (args.bc_list) {
        int val = 0;
        dict_read2(args.barcodes, args.bc_list, &val);
        if (dict_size(args.barcodes) == 0) error("Barcode list is empty?");
        args.fix_barcodes = 1;
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

        // detect output format from extension, default TSV
        const char *ext = strrchr(args.output_fname, '.');
        if (ext && strcmp(ext, ".csv") == 0) args.is_csv = 1;
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

    if (args.out && args.out != stdout) fclose(args.out);
}

int cov_main(int argc, char **argv)
{
    double t_real = realtime();

    if (parse_args(argc, argv)) return cov_usage();

    int n_cells = dict_size(args.barcodes);
    if (n_cells == 0) n_cells = 1024;
    int m_cells = n_cells;
    uint64_t *covered = calloc(m_cells, sizeof(uint64_t));
    if (covered == NULL) error("Failed to allocate memory.");

    int i;
    for (i = 0; i < args.hdr->n_targets; ++i) {
        int len = args.hdr->target_len[i];
        int last = 0;
        int end = last + ws1m;
        if (end > len) end = len;

        for (;;) {
            struct depth *d = bam2depth(args.idx, i, last, end, BED_STRAND_UNK,
                                        args.fp, args.mapq_thres, 0,
                                        args.barcodes, args.tag, args.umi_tag,
                                        1,             // split_by_tag = 1
                                        0, NULL,       // no alias
                                        args.fix_barcodes, -1);

            if (d) {
                struct depth *cur = d;
                while (cur) {
                    int id = cur->id;
                    if (id >= 0 && cur->dep1 + cur->dep2 >= args.min_depth) {
                        while (id >= m_cells) {
                            int old_m = m_cells;
                            m_cells = m_cells ? m_cells * 2 : 1024;
                            covered = realloc(covered, m_cells * sizeof(uint64_t));
                            memset(covered + old_m, 0, (m_cells - old_m) * sizeof(uint64_t));
                        }
                        covered[id]++;
                    }
                    cur = cur->next;
                }
                depth_destroy(d);
            }

            // handle newly discovered barcodes added by bam2depth
            int cur_n = dict_size(args.barcodes);
            while (cur_n > m_cells) {
                int old_m = m_cells;
                m_cells = m_cells ? m_cells * 2 : 1024;
                covered = realloc(covered, m_cells * sizeof(uint64_t));
                memset(covered + old_m, 0, (m_cells - old_m) * sizeof(uint64_t));
            }

            if (end == len) break;

            last = end;
            end = end + ws1m;
            if (end > len) end = len;
        }
    }

    int sep_char = args.is_csv ? ',' : '\t';
    write_output(covered, dict_size(args.barcodes), sep_char);

    free(covered);
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
