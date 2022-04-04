#include "utils.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "number.h"

static struct args {
    const char *input_fname;
    const char *output_fname;
    htsFile *in;
    FILE *out;
    int filter;
    int file_th;
    int fasta;
    int n_tag;
    char **tags;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .in = NULL,
    .out = NULL,
    .filter = 0,
    .file_th = 5,
    .fasta = 0,
    .n_tag = 0,
    .tags = NULL,
};

char **tag_split(const char *_s, int *n)
{
    kstring_t str = {0,0,0};
    kputs(_s, &str);
    int *s = ksplit(&str, ',', n);
    assert(*n>0);
    char **tags = malloc(sizeof(char*)*(*n));
    int j;
    for (j = 0; j < *n; ++j) {
        tags[j] = strdup(str.s+s[j]);
        if (strlen(tags[j]) != 2) error("Unknown tag format. Only two character allowed.");
    }
    free(str.s);
    free(s);
    return tags;
}
static int parse_args(int argc, char **argv)
{
    int i;
    const char *file_th = NULL;
    const char *tags = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if ( strcmp(a, "-h") == 0 || strcmp(a, "-help") == 0) return 1;
        else if (strcmp(a, "-f") == 0) {
            args.filter = 1;
            continue;
        }
        else if (strcmp(a, "-fa") == 0) {
            args.fasta = 1;
            continue;
        }
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-tags") == 0) var = &tags;
        
        if (var != 0) {
            if (argc == i) error("Miss an argument after %s.",a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;            
        }

        error("Unknown argument, %s", a);
    }

    if (args.input_fname == NULL) error("No input bam.");

    if (tags == NULL) error("No tag specified.");

    if (file_th) args.file_th = str2int((char*)file_th);

    args.tags = tag_split(tags, &args.n_tag);

    args.in = hts_open(args.input_fname, "r");
    CHECK_EMPTY(args.in, "Failed to open input bam.");

    htsFormat type = *hts_get_format(args.in);
    if (type.format != bam && type.format != sam)
        error("Unsupported input format, only support BAM/SAM/CRAM format.");

    hts_set_threads(args.in, args.file_th);

    args.out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    if (!args.out) error("%s : %s. ", args.output_fname, strerror(errno));
    
    return 0;
                                      
}
void memory_release()
{
    hts_close(args.in);
    fclose(args.out);
}

// reverse(), get_read(), and get_qual() copied from samtools/bam_fastq.c
// please refer to orginal resource if you use it
static char *reverse(char *str)
{
    int i = strlen(str)-1,j=0;
    char ch;
    while (i>j) {
        ch = str[i];
        str[i]= str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
}
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

/* return the read, reverse complemented if necessary */
static char *get_read(const bam1_t *rec)
{
    int len = rec->core.l_qseq + 1;
    char *read = calloc(1, len);
    char *seq = (char *)bam_get_seq(rec);
    int n;

    if (!read) return NULL;

    for (n=0; n < rec->core.l_qseq; n++) {
        if (rec->core.flag & BAM_FREVERSE) read[n] = seq_nt16_str[seq_comp_table[bam_seqi(seq,n)]];
        else                               read[n] = seq_nt16_str[bam_seqi(seq,n)];
    }
    if (rec->core.flag & BAM_FREVERSE) reverse(read);
    return read;
}

/*
 * get and decode the quality from a BAM record
 */
static int get_quality(const bam1_t *rec, char **qual_out)
{
    char *quality = calloc(1, rec->core.l_qseq + 1);
    char *q = (char *)bam_get_qual(rec);
    int n;

    if (!quality) return -1;

    if (*q == '\xff') {
        free(quality);
        *qual_out = NULL;
        return 0;
    }

    for (n=0; n < rec->core.l_qseq; n++) {
        quality[n] = q[n]+33;
    }
    if (rec->core.flag & BAM_FREVERSE) reverse(quality);
    *qual_out = quality;
    return 0;
}


extern int bam2fq_usage();
int bam2fq(int argc, char **argv)
{
    if (parse_args(argc, argv)) return bam2fq_usage();
    bam_hdr_t *hdr = sam_hdr_read(args.in);
    if (hdr == NULL) error("Failed to read header.");
    bam1_t *b = bam_init1();
    kstring_t str = {0,0,0};
    kstring_t name = {0,0,0};
    int ret;
    while ((ret = sam_read1(args.in, hdr, b)) >= 0) {
        name.l = 0;
        str.l = 0;
        int i;
        int filter = 0;
        for (i = 0; i < args.n_tag; ++i) {
            uint8_t *tag = bam_aux_get(b, args.tags[i]);
            if (!tag) {
                if (args.filter) {filter = 1; break;}
                continue;
            }
            else {
                kputs("|||", &name);
                kputs(args.tags[i], &name); kputc(':', &name);
                kputc(*(char*)tag, &name);  kputc(':', &name);
                kputs((char*)(tag+1), &name);
            }
        }
        
        if (filter == 1) continue;

        char *seq = get_read(b);
        
        if (args.fasta) {
            kputc('>', &str); kputs((char*)b->data, &str); kputs(name.s, &str); kputc('\n', &str);
            kputs(seq, &str);
            kputc('\n', &str);
        }
        else {
            char *qual = calloc(1, b->core.l_qseq+1);
            get_quality(b, &qual);

            kputc('@', &str); kputs((char*)b->data, &str); kputs(name.s, &str); kputc('\n', &str);
            kputs(seq, &str);
            // int i;
            //uint8_t *s = bam_get_seq(b);
            //for (i = 0; i < b->core.l_qseq; ++i) kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], &str);
            kputc('\n', &str);
            kputs("+\n", &str);
            //s = bam_get_qual(b);
            //if (s[0] == 0xff) for (i = 0; i < b->core.l_qseq; ++i) kputc('I', &str);
            //else for (i = 0; i < b->core.l_qseq; ++i) kputc(s[i] + 33, &str);
            kputs(qual, &str);            
            kputc('\n', &str);
            free(qual);
        }
        free(seq);
        fputs(str.s, args.out);
    }
    if (str.m) free(str.s);
    if (name.m) free(name.s);
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    memory_release();

    return 0;
}
