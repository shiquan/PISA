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
    int all_value;
    int mapq_thres;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .file_th = 4,
    .n_tag = 0,
    .print_rname = 0,
    .tags = NULL,
    .all_value = 0,
    .mapq_thres = 0,
};

extern int bam_extract_usage();

static int parse_args(int argc, char **argv)
{
    int i;
    const char *file_thread = NULL;
    const char *tags = NULL;
    const char *mapq = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        if (strcmp(a, "-tags") == 0) var = &tags;
        else if (strcmp(a, "-out") == 0 || strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-n") == 0) {
            args.print_rname = 1;
            continue;
        }
        else if (strcmp(a, "-all") == 0) {
            args.all_value = 1;
            continue;
        }
        else if (strcmp(a, "-q") == 0) var = &mapq;
        
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

    if (tags == 0) error("No tags specified.");
    if (file_thread) args.file_th = str2int((char*)file_thread);

    if (mapq) args.mapq_thres = str2int(mapq);
    
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
        if (b->core.flag & BAM_FQCFAIL) continue;
        if (b->core.flag & BAM_FSECONDARY) continue;
        if (b->core.qual < args.mapq_thres) continue;
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
            if (!tag) {
                if (args.all_value == 1) {
                    is_empty = 1;
                    break;
                }
                kputc('.', &str);
            }
            else {
                is_empty = 0;
                if (*tag == 'S' || *tag == 's' || *tag == 'c' || *tag == 'i' || *tag == 'I') {
                    int64_t va = bam_aux2i(tag);
                    kputw(va, &str);
                } else if (*tag == 'f' || *tag == 'd') {
                    double va = bam_aux2f(tag);
                    kputd(va, &str);
                } else if (*tag == 'H' || *tag == 'Z') {
                    char *va = bam_aux2Z(tag);
                    kputs(va, &str);
                }
                
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
