#include "utils.h"
#include "bed.h"
#include "gtf.h"

int gtf2bed_usage()
{
    fprintf(stderr, "# Convert GTF to BED.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m gtf2bed -o merged.bed in.gtf.gz\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  -o    [FILE]                          Output bed file.\n");
    fprintf(stderr, "  -type [gene|transcript|exon]          Covert to bed.\n");
    fprintf(stderr, "  -name [none|gene|transcript|exon]     Set name for bed.\n");
    fprintf(stderr, "\n");
    return 1;
}

static struct args {
    const char *input_fname;
    const char *output_fname;
    // 1 for gene, 2 for transcript, 3 for exon
    int type;

    // 0 for none, 1 for gene, 2 for transcript, 3 for exon
    int name_type;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .name_type = 0,
    .type = 0
};

static int parse_args(int argc, char **argv)
{
    if (argc == 1) return gtf2bed_usage();

    int i;
    const char *name_type = NULL;
    const char *type = NULL;
    
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a , "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-name") == 0) var = &name_type;
        else if (strcmp(a, "-type") == 0) var = &type;
        
        if (var != 0) {
            if (i == argc) error("missing an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (a[0] == '-' && a[1] != '\0') error("unknown argument, %s", a);
        
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("unknown argument, %s", a);
    }

    if (args.input_fname == NULL) error("No input specfied.");
    
    if (type == NULL) error("No type specified, shoule be gene/transcript/exon.");

    if (strcmp(type, "gene") == 0) args.type = 1;
    else if (strcmp(type, "transcript") == 0) args.type = 2;
    else if (strcmp(type, "exon") == 0) args.type = 3;
    else {
        error("unknown type, only support gene/transcript/exon");
    }
    
    if (name_type) {
        if (strcmp(name_type, "none") == 0) args.name_type = 0;
        else if (strcmp(name_type, "gene") == 0) args.name_type = 1;
        else if (strcmp(name_type, "transcript") == 0) args.name_type = 2;
        else if (strcmp(name_type, "exon") == 0) args.name_type = 3;
        else {
            error("unknown name type, only support none/gene/transcript/exon");
        }
    }
    return 0;
}
int gtf2bed_main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return gtf2bed_usage();

    struct gtf_spec *G = gtf_read(args.input_fname, 0x3);
    struct bed_spec *B = gtf2bed(G, NULL, args.type, args.name_type, 1);

    bed_spec_dedup(B, 1);
    bed_spec_write(B, args.output_fname, 0, 0);

    gtf_destroy(G);
    bed_spec_destroy(B);
    return 0;
}
