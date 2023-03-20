#include "utils.h"
#include "bed.h"

static struct args {
    int n_file;
    const char **input_fnames;
    const char *output_fname;
    int upstream;
    int downstream;
    int stranded;
    int check_name;
} args = {
    .n_file = 0,
    .input_fnames = NULL,
    .output_fname = NULL,
    .upstream = 0,
    .downstream = 0,
    .stranded = 1,
    .check_name = 0
};

static int mergebed_usage()
{
    fprintf(stderr, "# Merge overlaped regions in bed files.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m mergebed -o merged.bed sample1.bed sample2.bed\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m mergebed -up 500 -down 500 -o flank.bed peaks.bed\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-o    [FILE]    Output bed file.\n");
    fprintf(stderr, "-s              Ignore strand.\n");
    fprintf(stderr, "-up   [INT]     Enlarge regions upstream.\n");
    fprintf(stderr, "-down [INT]     Enlarge regions downstream.\n");
    fprintf(stderr, "-name           Merge regions by bed name.\n");
    
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * This tool accepts 3 columns or 6 columns bed file(s), strand (+/-) is set in column 6.\n");
    fprintf(stderr, " * By default, merging is done with respect to strandness unless -s is set.\n");
    fprintf(stderr, " * -up/-down is set respect to strandness, so upstream of plus strand is downstream of minus strand.\n");
    
    fprintf(stderr, "\n");
    return 1;
}

static int parse_args(int argc, char **argv)
{
    int i;
    const char *upstream = 0;
    const char *downstream = 0;
    
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-up") == 0) var = &upstream;
        else if (strcmp(a, "-down") == 0) var = &downstream;
        else if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-s") == 0) {
            args.stranded = 0;
            continue;
        }
        else if (strcmp(a, "-name") == 0) {
            args.check_name = 1;
            continue;
        }
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1] != '\0') error("Unknown argument, %s",a);

        if (args.input_fnames == NULL) {
            args.n_file = 1;
            args.input_fnames = malloc(sizeof(char*));
            args.input_fnames[0] = a;
        } else {
            args.input_fnames = realloc(args.input_fnames, (args.n_file+1)*sizeof(char*));
            args.input_fnames[args.n_file] = a;
            args.n_file++;
        }        
    }

    if (args.n_file == 0) error("No input bed file(s).");
    
    return 0;
}

static void memory_release()
{
    free(args.input_fnames);
}

int mergebed(int argc, char **argv)
{
    if (parse_args(argc, argv)) return mergebed_usage();

    int i;
    struct bed_spec *B = bed_spec_init();
    
    for (i = 0; i < args.n_file; ++i) {
        B = bed_read0(B, args.input_fnames[i]);
    }

    bed_spec_merge1(B, args.stranded, args.upstream, args.downstream, 0, args.check_name);
    
    bed_spec_write(B, args.output_fname, 0);

    bed_spec_destroy(B);

    memory_release();
    
    return 0;
}
