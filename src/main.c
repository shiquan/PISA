#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int usage()
{
    fprintf(stderr, "SingleCellTools - Collection of tools to process single cell omics data.\n");
    fprintf(stderr, "--- FASTQ ---\n");
    fprintf(stderr, " parse      Parse cell barcode, sample barcode and UMI from fastq reads to read name.\n");
    fprintf(stderr, " trim       Trim TN5 mosic ends, poly As.\n");
    fprintf(stderr, "--- BAM ---\n");
    //fprintf(stderr, " sam2bam\n");
    //fprintf(stderr, " rmdup\n");
    fprintf(stderr, " anno       Annotate peak or gene names into BAM attributions.\n");
    fprintf(stderr, " attrcnt    Count raw reads and reads with predefined tag for eacg cell barcode.\n");
    //fprintf(stderr, " count\n");
    fprintf(stderr, "\n");
    return 1;
}
int main(int argc, char *argv[])
{
    extern int bam_anno_attr(int argc, char *argv[]);
    extern int bam_count_attr(int argc, char *argv[]);
    if (argc == 1) return usage();
    else if ( strcmp(argv[1], "anno") == 0 ) return bam_anno_attr(argc-1, argv+1);
    else if ( strcmp(argv[1], "attrcnt") == 0) return bam_count_attr(argc-1, argv+1);
    else return usage();
    return 0;
}
