#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "number.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *str;
    int mapq_thres;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .str = NULL,
    .mapq_thres = 0,
};

static int usage()
{
    fprintf(stderr, "# Just add tags to reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m addtag -str CB:Z:CELL1,LB:Z:PoII2 -q 20 -o out.bam in.bam\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-o    [FILE]    Output bam file.\n");
    fprintf(stderr, "-str  [string]  TAGs.\n");
    fprintf(stderr, "-@    [INT]     Threads to pack file.\n");
    fprintf(stderr, "-mapq [INT]     Mapping quality score to filter reads.\n");
    fprintf(stderr, "\n");
    return 1;
}
int add_tags(int argc, char **argv)
{
    if (argc == 1) return usage();
    
    int i;
    const char *mapq = NULL;
    const char *threads = NULL;
    int n_thread = 1;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0) return usage();
        
        if (strcmp(a, "-str") == 0) var = &args.str;
        else if (strcmp(a, "-o") == 0 || strcmp(a, "-out") == 0)
            var = &args.output_fname;
        else if (strcmp(a, "-mapq") == 0 || strcmp(a, "-q") == 0)
            var = &mapq;
        else if (strcmp(a, "-@") == 0)
            var = &threads;
        
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

    if (args.str == NULL) error("No -str specified.");
    if (mapq) args.mapq_thres = str2int(mapq);
    if (threads) n_thread = str2int(threads);

    if (args.output_fname == NULL) error("No output.");
    if (args.input_fname  == NULL) error("No input.");
    
    kstring_t str = {0,0,0};

    kputs(args.str, &str);

    int n_tag = 0;
    int *s = ksplit(&str, ',', &n_tag);

    for (i = 0; i < n_tag; ++i) {
        char *tag = str.s +s[i];
        if (strlen(tag) < 6) error("Unknown tag format %s", tag);
        if (tag[2] != ':' && tag[4] != ':') error("Unknown tag format %s", tag);
        tag[2] = 0;
        tag[4] = 0;
    }
    htsFile *in = hts_open(args.input_fname, "r");

    if (in == NULL) error("%s : %s.", args.input_fname, strerror(errno));

    bam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == NULL) error("Failed to read header.");

    
    htsFile *out = hts_open(args.output_fname, "wb");

    hts_set_threads(out, n_thread);
    
    if (sam_hdr_write(out, hdr)) error("Failed to write SAM header.");
    
    int ret;
    bam1_t *b;
    
    b = bam_init1();

    while ((ret = sam_read1(in, hdr, b)) >= 0) {

        if (b->core.qual < args.mapq_thres) continue;

        int i;
        for (i = 0; i < n_tag; ++i) {
            uint8_t *data;
            char *tag = str.s+s[i];
            
            if ((data = bam_aux_get(b, tag)) != NULL) bam_aux_del(b, data);
            //debug_print("%s\t%c\t%d\t%s", tag, tag[3], strlen(tag+5), tag+5);
            bam_aux_append(b, tag, tag[3], strlen(tag+5)+1, (uint8_t*)(tag+5));
        }
        
        if (sam_write1(out, hdr, b) == -1)
            error("Failed to write SAM.");
    }

    if (str.m) free(str.s);
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);
    sam_close(out);
    free(s);
    return 0;
}
   
