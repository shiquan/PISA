// Count all reads and reads with predefined tag for each cell barcode in the bam file
#include "utils.h"
#include "number.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/kseq.h"
#include "bed_lite.h"
#include "barcode_list.h"
#include <zlib.h>

static int usage()
{
    fprintf(stderr, "BamCount -cb CR -tag PK,TS -o count.txt -list barcode.txt in.bam\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *barcode_fname;
    const char *cb_tag; // attribute in BAM
    int n_tag;
    char **tags;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .barcode_fname = NULL,
    .cb_tag = NULL,
    .n_tag = 0,
    .tags = NULL,
};
static int parse_args(int argc, char **argv)
{
    int i;
    const char *tag = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-cb") == 0) var = &args.cb_tag;
        else if (strcmp(a, "-tag") == 0) var = &tag;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
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
    CHECK_EMPTY(args.output_fname, "-o must be set.");
    CHECK_EMPTY(args.cb_tag, "-cb must be set.");
    CHECK_EMPTY(args.input_fname, "Input bam must be set.");
    if (tag) {
        kstring_t str = {0,0,0};
        kputs(tag, &str);
        int n = 0;
        int *s = ksplit(&str, ',', &n);
        args.tags = malloc(n*sizeof(char*));
        args.n_tag = n;
        int j;
        for (j = 0; j < n; ++j) {
            args.tags[j] = strdup(str.s + s[j]);
        }
        free(s);
        free(str.s);
    }
    return 0;
}
struct count {
    uint32_t raw;
    uint32_t *c; // counts for tags
};
static struct count *count_init(int n)
{
    struct count *c = malloc(sizeof(*c));
    //memset(c, 0, sizeof(*c));
    c->c = calloc(n, sizeof(uint32_t));
    c->raw = 0;
    return c;
}
int bam_count_attr(int argc, char *argv[])
{
    double t_real;
    t_real = realtime();
    if (parse_args(argc, argv)) return usage();
    
    htsFile *fp  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    CHECK_EMPTY(hdr, "Failed to open header.");
    FILE *out = fopen(args.output_fname, "w");
    CHECK_EMPTY(out, "%s : %s.", args.output_fname, strerror(errno));
    fprintf(out, "CELL_BARCODE\tReadCounts");
    int i;
    for (i = 0; i < args.n_tag; ++i) fprintf(out, "\t%s", args.tags[i]);
    fputc('\n', out);
    int is_dyn_alloc = 1;

    struct lbarcode *barcode = barcode_init();
    if (args.barcode_fname) {
        barcode_read(barcode, args.barcode_fname);
        is_dyn_alloc = 0;
    }

    bam1_t *b;
    bam1_core_t *c;
    int ret;
    b = bam_init1();
    c = &b->core;

    while ((ret = sam_read1(fp, hdr, b)) >= 0) {
        //char *name = hdr->target_name[c->tid];
        uint8_t *tag = bam_aux_get(b, args.cb_tag);
        if (!tag) error("No cell barcode tag at alignment. %s:%d", hdr->target_name[c->tid], c->pos+1);
        char *name = (char*)(tag+1);
        int id = -1;
        if (is_dyn_alloc == 0) {
            id = barcode_select(barcode, name);
            if (id == -1) continue;
        }
        else {
            id = barcode_push(barcode, name);
        }
        struct barcode *d = &barcode->b[id];
        if (d->data == NULL) {
            struct count *t = count_init(args.n_tag);
            d->data = (void*)t;
        }

        struct count *data = (struct count*)d->data;
        data->raw++;
        for (i = 0; i < args.n_tag; ++i) {
            uint8_t *g = bam_aux_get(b, args.tags[i]);
            if (g) data->c[i]++;
        }
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    if (ret != -1) warnings("Truncated file?");
    
    int j;
    for (i = 0; i < barcode->n; ++i) {
        struct barcode *d = &barcode->b[i];
        if (d->data == NULL) continue; // if not exists, skip
        struct count *data = (struct count*)d->data;
        fprintf(out, "%s\t%u", d->s, data->raw);
        for (j = 0; j < args.n_tag; ++j)
            fprintf(out, "\t%u", data->c[j]);
        fputc('\n', out);
        if (data->c) free(data->c);
        free(data);        
    }
    barcode_destory(barcode);
    fclose(out);
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;    
}


        
