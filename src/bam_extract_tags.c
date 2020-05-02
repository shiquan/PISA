#include "utils.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "number.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    int file_th;
    int n_tag;
    int print_rname;
    char **tags;    
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .file_th = 4,
    .n_tag = 0,
    .print_rname = 0,
    .tags = NULL,
};

extern int bam_extract_usage();

static int parse_args(int argc, char **argv)
{
    int i;
    const char *file_thread = NULL;
    const char *tags = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tags") == 0) var = &tags;        
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-n") == 0) {
            args.print_rname = 1;
            continue;
        }
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
    // if (args.output_fname == 0) error("No output file specified.");
    if (tags == 0) error("No tags specified.");
    if (file_thread) args.file_th = str2int((char*)file_thread);

    kstring_t str = {0,0,0};
    kputs(tags, &str);
    int *s = ksplit(&str, ',', &args.n_tag);
    assert(args.n_tag>0);
    args.tags = malloc(sizeof(char*)*args.n_tag);
    int j;
    for (j = 0; j < args.n_tag; ++j) {
        args.tags[j] = strdup(str.s+s[j]);
        if (strlen(args.tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
    }
    free(str.s);
    free(s);
    
    return 0;
}

int bam_extract_tags(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bam_extract_usage();
    htsFile *in  = hts_open(args.input_fname, "r");
    CHECK_EMPTY(in, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    bam_hdr_t *hdr = sam_hdr_read(in);
    CHECK_EMPTY(hdr, "Failed to open header.");
    hts_set_threads(in, args.file_th);

    FILE *out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    if (out == NULL) error("%s : %s.", args.output_fname, strerror(errno));

    int ret;
    bam1_t *b;
    b = bam_init1();
    kstring_t str = {0,0,0};
    int is_empty;
    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        int i;
        is_empty = 1;
        str.l = 0;
        if (args.print_rname) {
            kputs((char*)b->data, &str);
            kputc('\t', &str);
        }
        for (i = 0; i < args.n_tag; ++i) {
            if (i) kputc('\t', &str);
            uint8_t *tag = bam_aux_get(b, args.tags[i]);
            if (!tag) kputc('.', &str);
            else {
                is_empty = 0;
                if (*tag == 'A') kputc(tag[1], &str);
                else kputs((char*)(tag+1), &str);
            }
        }
        kputc('\n', &str);
        if (is_empty==0) fputs(str.s, out);
    }
    if (str.m) free(str.s);
    bam_destroy1(b);
    
    fclose(out);
    bam_hdr_destroy(hdr);
    sam_close(in);

    int i;
    for (i = 0; i < args.n_tag; ++i) free(args.tags[i]);
    free(args.tags);
    return 0;
}
