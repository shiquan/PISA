#include "utils.h"
#include "gtf.h"
#include "htslib/kstring.h"

int gtf_format_usage()
{
    fprintf(stderr, "# format and order GTF file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m gtffmt in.gtf\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, " -o    [FILE]    Output GTF file\n");
    fprintf(stderr, " -f              Only export gene, transcript, and exon records.\n");
    fprintf(stderr, " -key  [all]     Export keys in optional fields, export all keys in default.\n");
    fprintf(stderr, "\n");
    return 1;
}

int gtf_format(int argc, char **argv)
{
    if (argc == 1) return gtf_format_usage();
    int only_exon = 0;
    int i;
    const char *input_fname = NULL;
    const char *output_fname = NULL;
    const char *key_str = NULL;
    
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0) return gtf_format_usage();
        
        if (strcmp(a, "-o") == 0 || strcmp(a, "-out") == 0)
            var = &output_fname;
        else if(strcmp(a, "-f") == 0) {
            only_exon = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        
        if (input_fname == NULL) {
            input_fname = a;
            continue;
        }
        error("Unknown argument, %s", a);
    }

    if (output_fname == NULL) error("No output.");
    if (input_fname  == NULL) error("No input.");

    struct dict *keys = NULL;
    if (key_str) {
        int *s;
        int n;
        kstring_t str = {0,0,0};
        kputs(key_str, &str);
        s = ksplit(&str,',', &n);
        int i;
        keys = dict_init();
        for (i = 0; i < n; ++i) {
            char *key = str.s + s[i];
            dict_push(keys, key);
        }
    }
    struct gtf_spec *G = gtf_read(input_fname, only_exon);
    gtf_dump(G, output_fname, keys);

    gtf_destroy(G);
    if (keys) dict_destroy(keys);
    return 0;
}
