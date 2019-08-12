#include "utils.h"
#include "fastq.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "number.h"

KSEQ_INIT(gzFile, gzread)

static int usage()
{
    fprintf(stderr, "fastq_split [options] in.fq\n");
    fprintf(stderr, "      -chunk     Minimize records per file.\n");
    fprintf(stderr, "      -out       Outdir.\n");
    fprintf(stderr, "      -tag       Tags.\n");
    fprintf(stderr, "\n");
    return 1;
}
static struct args {
    const char *outdir;
    const char *input_fname;
    int n_tag;
    char **tags;
    int chunk;    
} args = {
    .outdir = "./",
    .input_fname = NULL,
    .n_tag = 0,
    .tags = NULL,
    .chunk = 1000,
};
static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *tags = NULL;
    const char *chunk = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-out") == 0 || strcmp(a, "-o") == 0) var = &args.outdir;
        else if (strcmp(a, "-chunk") == 0) var = &chunk;
        if (var != 0) {
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        
        error("Unknown argument : %s", a);
    }
    // todo: check file, NO streaming allowed
    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (tags == NULL) error("-tag must be set.");

    kstring_t str = {0,0,0};
    kputs(tags, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    args.n_tag = n;
    args.tags = malloc(n*sizeof(char*));
    for (i = 0; i <n; ++i) args.tags[i] = strdup(str.s+s[i]);
    free(s); free(str.s);
    if (chunk) args.chunk = str2int((char*)chunk);
    if (args.chunk < 0) args.chunk = 1000;
    return 0;
}

int LFR_split(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();
    gzFile fp = gzopen(args.input_fname, "r");
    if (fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    kseq_t *ks;
    ks = kseq_init(fp);
    if (ks == NULL) error("Failed to load fastq.");
    char *last_name = NULL;
    kstring_t file = {0,0,0};
    int export_last = 0;
    FILE *out = NULL;
    int idx = 0;
    int file_idx = 0;
    ksprintf(&file, "%s/split_%d.fq", args.outdir, file_idx);
    out = fopen(file.s, "w");
    if (out == NULL) error("%s : %s", file.s, strerror(errno));
    for (;;) {
        if (export_last) {            
            export_last = 0;
            goto export_record;
        }
        if (kseq_read(ks) < 0) break;
        idx++;
        if (idx >= args.chunk) {
            char **names = fastq_name_pick_tags(ks->name.s, args.n_tag, args.tags);
            kstring_t str = {0,0,0};
            int i;
            for (i = 0; i < args.n_tag; ++i) {
                if (names == NULL || names[i] == NULL) error("No tag at read name. %s", ks->name.s);
                kputs(names[i], &str);
                free(names[i]);
            }
            free(names);
            if (last_name) {
                if (strcmp(last_name, str.s) == 0) {
                    free(str.s);
                    goto export_record;
                }
                else {
                    free(last_name); last_name = NULL;
                    free(str.s);
                    last_name = NULL;
                    export_last = 1;
                    fclose(out);
                    file_idx++;
                    idx = 0;
                    file.l = 0;
                    ksprintf(&file, "%s/split_%d.fq", args.outdir, file_idx);
                    out = fopen(file.s, "w");
                    if (out == NULL) error("%s : %s", file.s, strerror(errno));
                    continue;
                }
            }
            else last_name = str.s;
        }
        
      export_record:
        fprintf(out, "@%s\n%s\n+\n%s\n", ks->name.s, ks->seq.s, ks->qual.s);        
    };
    fclose(out);
    if (last_name) free(last_name);
    return 0;
}

int main(int argc, char **argv)
{
    return LFR_split(argc, argv);
}
