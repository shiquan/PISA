#include "utils.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv)
{
    int length = 100;
    const char *in = NULL;
    const char *len_str = NULL;
    int i;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0) {
            fprintf(stderr, "seqcut -length 100 in.fq\n");
            return 1;
        }

        if (strcmp(a, "-length") == 0) var = &len_str;


        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (in == 0) {
            in = a;
            continue;
        }
        error("Unknown argument: %s", a);
    }

    if (in == NULL) error("No input.");
    if (len_str) length = atoi(len_str);
    if (length < 0) length = 100;
    gzFile fp = gzopen(in, "r");
    if (fp == NULL) error("%s : %s.", in, strerror(errno));
    kseq_t *ks = kseq_init(fp);
    kstring_t str ={0,0,0};
    int ret;
    while ((ret = kseq_read(ks)) >= 0) {
        if (ks->seq.l >= length) {
            str.l = 0;
            int i;
            for (i = 0; i + length <= ks->seq.l; i+=50) {
                if (ks->qual.l) kputc('@', &str);
                else kputc('>',&str);
                kputs(ks->name.s, &str);
                kputc('\n', &str);
                kputsn(ks->seq.s+i, length, &str);
                kputc('\n', &str);
                if (ks->qual.l) {
                    kputs("+\n", &str);
                    kputsn(ks->qual.s+i, length, &str);
                    kputc('\n', &str);
                }
            }
            puts(str.s);
        }
    }
    if (str.m) free(str.s);
    kseq_destroy(ks);
    gzclose(fp);
    return 0;
}
