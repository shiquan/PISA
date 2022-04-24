#include "utils.h"
#include "version.h"
#include "pisa_version.h"
#include <string.h>
int usage()
{
    fprintf(stderr, "\n\x1b[1mPISA\x1b[0m - a collection of tools for single cell data pre-processing and interpretation.\n");
    fprintf(stderr, "Version: %s + htslib: %s\n", PISA_VERSION, HTS_VERSION_TEXT);
    fprintf(stderr, "Contact: Quan SHI [shiquan(AT)genomics.cn]\n");
    fprintf(stderr, "\nCommands:\n");
    fprintf(stderr, "\n--- Processing FASTQ/FASTQ+\n");
    fprintf(stderr, "    parse2     A new generic tool to parse barcodes from FASTQ.\n");
    fprintf(stderr, "    parse      Parse barcodes from FASTQ reads to FASTQ+.\n");
    //fprintf(stderr, "    pick       Pick FASTQ+ records with tags.\n");
    //fprintf(stderr, "    attrcnt    Count raw reads and tag values per cell for FASTQ+.\n");
    fprintf(stderr, "    fsort      Sort FASTQ+ records by barcodes.\n");
    fprintf(stderr, "    stream     Perform user-defined process for each read block.\n");
    fprintf(stderr, "    addtags    Add tag string to FASTQ reads.\n");
    
    fprintf(stderr, "\n--- Processing BAM\n");
    fprintf(stderr, "    sam2bam    Parse FASTQ+ read name and convert SAM to BAM.\n");
    fprintf(stderr, "    rmdup      Remove PCR duplicates per molecular.\n");
    fprintf(stderr, "    pick       Pick alignments with tags.\n");
    fprintf(stderr, "    anno       Annotate functional regions or gene names.\n");
    fprintf(stderr, "    corr       Correct error prone UMIs. 1 mismatch considered.\n");
    fprintf(stderr, "    attrcnt    Count raw reads and tag values per cell.\n");
    fprintf(stderr, "    extract    Extract tag value from BAM.\n");
    fprintf(stderr, "    count      Count feature X cell matrix from BAMs.\n");
    fprintf(stderr, "    bam2fq     Convert BAM to FASTQ+ file with selected tags.\n");
    fprintf(stderr, "    bam2frag   Generate fragment file.\n");
    fprintf(stderr, "    depth      Coverage depth/UMI for target regions.\n");
    fprintf(stderr, "    addtags    Add tag string to BAM alignments.\n");

    fprintf(stderr, "\n--- Processing fragment file.\n");
    fprintf(stderr, "    count2     Count peak X cell matrix from fragment file.\n");

    fprintf(stderr, "\n--- Some experimental ideas. Not stable, just for test.\n");
    fprintf(stderr, "    fusion     Predict gene fusion based on UMIs. **experiment**\n");
    fprintf(stderr, "    gtffmt     Format and reorder GTF file. **experiment**\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Use `\x1b[1mPISA\x1b[0m command -h` for help information.\n");
    fprintf(stderr, "\n");
    return 1;
}
int main(int argc, char *argv[])
{
    // process FQ
    extern int fastq_parse_barcodes(int argc, char *argv[]);
    extern int fastq_parse2(int argc, char **argv);
    extern int fsort(int argc, char **argv);
    extern int fastq_stream(int argc, char **argv);
    
    // process BAM
    extern int sam2bam(int argc, char *argv[]);
    extern int bam_rmdup(int argc, char *argv[]);
    extern int bam_anno_attr(int argc, char *argv[]);
    extern int bam_count_attr(int argc, char *argv[]);
    extern int bam_pick(int argc, char *argv[]);
    extern int bam_corr_umi(int argc, char **argv);
    extern int bam_extract_tags(int argc, char **argv);
    extern int count_matrix(int argc, char *argv[]);
    extern int bam2fq(int argc, char *argv[]);
    extern int bam2frag(int argc, char **argv);
    extern int gene_fusion(int argc, char **argv);
    extern int depth_main(int argc, char **argv);
    extern int add_tags(int argc, char **argv);

    // process fragment
    extern int fragment_count(int argc, char **argv);

    // process GTF
    extern int gtf_format(int argc, char **argv);
    
    if (argc == 1) return usage();
    else if (strcmp(argv[1], "parse") == 0) return fastq_parse_barcodes(argc-1, argv+1);
    else if (strcmp(argv[1], "parse2") == 0) return fastq_parse2(argc-1, argv+1);
    else if (strcmp(argv[1], "fsort") == 0) return fsort(argc-1, argv+1);
    else if (strcmp(argv[1], "stream") == 0) return fastq_stream(argc-1, argv+1);
    else if (strcmp(argv[1], "sam2bam") == 0) return sam2bam(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2fq") == 0) return bam2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "rmdup") == 0) return bam_rmdup(argc-1, argv+1);
    else if (strcmp(argv[1], "anno") == 0) return bam_anno_attr(argc-1, argv+1);
    else if (strcmp(argv[1], "corr") == 0) return bam_corr_umi(argc-1, argv+1);
    else if (strcmp(argv[1], "attrcnt") == 0) return bam_count_attr(argc-1, argv+1);
    else if (strcmp(argv[1], "extract") == 0) return bam_extract_tags(argc-1, argv+1);
    else if (strcmp(argv[1], "pick") == 0) return bam_pick(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2frag") == 0) return bam2frag(argc-1, argv+1);
    else if (strcmp(argv[1], "count") == 0) return count_matrix(argc-1, argv+1);
    else if (strcmp(argv[1], "count2") == 0) return fragment_count(argc-1, argv+1);
    else if (strcmp(argv[1], "fusion") == 0) return gene_fusion(argc-1, argv+1);
    // else if (strcmp(argv[1], "genecov") == 0) return gene_cov(argc-1, argv+1);
    // else if (strcmp(argv[1], "assem") == 0)  return fastq_assem(argc-1, argv+1);
    // else if (strcmp(argv[1], "segment") == 0) return fastq_segment(argc-1, argv+1);
    // else if (strcmp(argv[1], "segment2") == 0) return check_segment2(argc-1, argv+1);
    // else if (strcmp(argv[1], "cleanup") == 0) return LFR_cleanup(argc-1, argv+1);
    // else if (strcmp(argv[1], "overlap") == 0) return fastq_overlap(argc-1, argv+1);
    // else if (strcmp(argv[1], "impute") == 0) return LFR_impute(argc-1, argv+1);
    else if (strcmp(argv[1], "depth") == 0) return depth_main(argc-1, argv+1);
    else if (strcmp(argv[1], "addtags") == 0) return add_tags(argc-1, argv+1);
    else if (strcmp(argv[1], "gtffmt") == 0) return gtf_format(argc-1, argv+1);
    else return usage();
    return 0;
}
