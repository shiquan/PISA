#include "utils.h"
#include "bed.h"
#include "gtf.h"


int bed_flatten_usage()
{
    fprintf(stderr, "# Convert overlapped bed to flattern bed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m flatten -o merged.bed overlapped.bed\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-o    [FILE]    Output bed file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "For example:\n");
    fprintf(stderr, "reg1 ===========\n");
    fprintf(stderr, "reg2       ===========\n");
    fprintf(stderr, "flattening of regions:\n");
    fprintf(stderr, "reg1 ======\n");
    fprintf(stderr, "reg2       =====\n");
    fprintf(stderr, "reg3            ======\n");
    return 1;
}
static struct args {
    const char *input_fname;
    const char *output_fname;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
};

static int parse_args(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a , "-o") == 0) var = &args.output_fname;

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

    if (args.input_fname == NULL) return 1;   
    return 0;
}
int bed_flatten(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bed_flatten_usage();
    struct bed_spec *B = bed_read(args.input_fname);
    struct bed_spec *flatten = bed_spec_flatten(B);

    bed_spec_write(flatten, args.output_fname, 0, 0);
    bed_spec_destroy(B);
    bed_spec_destroy(flatten);

    return 0;
}
