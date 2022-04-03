#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "number.h"
#include "htslib/kseq.h"


static struct args {
    const char *input_fname;
    const char *output_fname;
    int mapq_thres;
    htsFile *in;
    int n_tag;
    kstring_t *tag_str;
    int *s;
    int n_thread;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .mapq_thres = 0,
    .in = NULL,
    .n_tag = 0,
    .tag_str = NULL,
    .s = NULL,
    .n_thread = 1,
};

static int usage()
{
    fprintf(stderr, "# Just add tags to reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m addtags -str CB:Z:CELL1,LB:Z:PoII2 -o out.bam in.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m addtags -str CB:Z:CELL1,LB:Z:PoII2 -o out.fastq in.fastq\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-o    [FILE]    Output file, bam or fastq, depends on input format\n");
    fprintf(stderr, "-str  [string]  TAGs.\n");
    fprintf(stderr, "-@    [INT]     Threads to pack file.\n");
    fprintf(stderr, "-mapq [INT]     Mapping quality score to filter mapped reads.\n");
    fprintf(stderr, "\n");
    return 1;
}


void check_tagstr(const char *tag_str)
{
    args.tag_str = malloc(sizeof(kstring_t));
    
    args.tag_str->s = NULL;
    args.tag_str->m = args.tag_str->l = 0;
    
    kputs(tag_str, args.tag_str);
    
    args.s = ksplit(args.tag_str, ',', &args.n_tag);

    if (args.n_tag < 1) error("Unknown tag format %s", tag_str);
    int i;
    for (i = 0; i < args.n_tag; ++i) {
        char *tag = args.tag_str->s + args.s[i];
        if (strlen(tag) < 6) error("Unknown tag format %s", tag);
        if (tag[2] != ':' && tag[4] != ':') error("Unknown tag format %s", tag);
        tag[2] = 0;
        tag[4] = 0;
    }
}
void generate_bamout()
{
    bam_hdr_t *hdr = sam_hdr_read(args.in);
    if (hdr == NULL) error("Failed to read header.");
    
    htsFile *out = hts_open(args.output_fname, "wb");

    hts_set_threads(out, args.n_thread);
    
    if (sam_hdr_write(out, hdr)) error("Failed to write SAM header.");
    
    int ret;
    bam1_t *b;
    
    b = bam_init1();
    
    while ((ret = sam_read1(args.in, hdr, b)) >= 0) {
        
        if (b->core.qual < args.mapq_thres) continue;

        int i;
        for (i = 0; i < args.n_tag; ++i) {
            uint8_t *data;
            char *tag = args.tag_str->s + args.s[i];
            
            if ((data = bam_aux_get(b, tag)) != NULL) bam_aux_del(b, data);
            bam_aux_append(b, tag, tag[3], strlen(tag+5)+1, (uint8_t*)(tag+5));
        }
        
        if (sam_write1(out, hdr, b) == -1)
            error("Failed to write SAM.");
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(out);    
}

#include <zlib.h>
#include "htslib/bgzf.h"

KSEQ_INIT(gzFile, gzread)
    ;
void generate_fqout()
{
    gzFile r = gzopen(args.input_fname, "r");
    BGZF *o = bgzf_open(args.output_fname, "w");

    bgzf_mt(o, args.n_thread, 256);
        
    kseq_t *ks = kseq_init(r);

    kstring_t str = {0,0,0};

    kstring_t tag_str = {0,0,0};

    // todo: use fname_update_tags() instead
    int i;
    for (i = 0; i < args.n_tag; ++i) {
        char *tag = args.tag_str->s + args.s[i];
        tag[2] = ':';
        tag[4] = ':';
        kputs("|||", &tag_str);
        kputs(tag, &tag_str);    
    }
    
    int ret;
    
    while( (ret = kseq_read(ks)) >= 0) {
        kputc('@', &str);
        kputs(ks->name.s, &str);
        kputs(tag_str.s, &str);
        kputc('\n', &str);
        kputs(ks->seq.s, &str);
        kputs("\n+\n", &str);
        kputs(ks->qual.s, &str);
        kputc('\n', &str);

        if (bgzf_write(o, str.s, str.l) < 0) error("Failed to write file.");
        str.l = 0;
    }

    free(str.s);
    free(tag_str.s);
    kseq_destroy(ks);
    gzclose(r);
    bgzf_close(o);
}

int add_tags(int argc, char **argv)
{
    if (argc == 1) return usage();
    
    int i;
    const char *mapq = NULL;
    const char *threads = NULL;
    const char *str = NULL;
    
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0) return usage();
        
        if (strcmp(a, "-str") == 0) var = &str;
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

    if (str == NULL) error("No -str specified.");
    if (mapq) args.mapq_thres = str2int(mapq);
    if (threads) args.n_thread = str2int(threads);

    if (args.output_fname == NULL) error("No output.");
    if (args.input_fname  == NULL) error("No input.");

    check_tagstr(str);
    
    args.in = hts_open(args.input_fname, "r");
    if (args.in == NULL) error("%s : %s.", args.input_fname, strerror(errno));

    htsFormat type = *hts_get_format(args.in);

    if (type.format != bam && type.format != fastq_format)
        error("Unsupported input format, only support BAM/FASTQ format.");
    
    if (type.format == bam) generate_bamout();
    else generate_fqout(); 
    
    if (args.tag_str->m) free(args.tag_str->s);
    free(args.tag_str);
    
    hts_close(args.in);
    free(args.s);
    
    return 0;
}
