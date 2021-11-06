#include "utils.h"
#include "pisa_version.h"
#include "version.h"
#include <string.h>

int usage()
{
    fprintf(stderr, "\nPISA - a preprocessing and interative suite for single-cell data analysis\n");
    fprintf(stderr, "Version: %s + htslib: %s\n", PISA_VERSION, HTS_VERSION_TEXT);
    fprintf(stderr, "Contact: Quan SHI [shiquan(AT)genomics.cn]\n");
    fprintf(stderr, "\nCommands:\n");
    fprintf(stderr, "\n--- Processing FASTQ\n");
    fprintf(stderr, "    parse      Parse barcodes from fastq reads.\n");
    fprintf(stderr, "    fsort      Sort fastq records by barcodes. **experiment**\n");
    
    fprintf(stderr, "\n--- Processing BAM\n");
    fprintf(stderr, "    sam2bam    Parse FASTQ+ read name and convert SAM to BAM.\n");
    fprintf(stderr, "    rmdup      Remove PCR duplicates per molecular.\n");
    fprintf(stderr, "    pick       Pick alignments with tags.\n");
    fprintf(stderr, "    anno       Annotate functional regions or gene names.\n");
    fprintf(stderr, "    corr       Correct error prone UMIs. 1 mismatch considered.\n");
    fprintf(stderr, "    attrcnt    Count raw reads and tag values per cell.\n");
    fprintf(stderr, "    extract    Extract tag value from BAM.\n");
    fprintf(stderr, "    count      Count matrix.\n");
    fprintf(stderr, "    bam2fq     Convert BAM to FASTQ+ file with selected tags.\n");
    fprintf(stderr, "    bam2frag   Generate fragment file.\n");
    fprintf(stderr, "    fusion     Predict gene fusion based on UMIs. **experiment**\n");
    fprintf(stderr, "\n");
    return 1;
}
int main(int argc, char *argv[])
{
    // process FQ
    extern int fastq_prase_barcodes(int argc, char *argv[]);
    //extern int fastq_trim_adaptors(int argc, char *argv[]);
    extern int fsort(int argc, char ** argv);

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
    // extern int gene_cov(int argc, char **argv);
    extern int bam2frag(int argc, char **argv);
    extern int gene_fusion(int argc, char **argv);


    if (argc == 1) return usage();
    else if (strcmp(argv[1], "parse") == 0) return fastq_prase_barcodes(argc-1, argv+1);
    //else if (strcmp(argv[1], "trim") == 0) return fastq_trim_adaptors(argc-1, argv+1);
    else if (strcmp(argv[1], "fsort") == 0) return fsort(argc-1, argv+1);
    else if (strcmp(argv[1], "sam2bam") == 0) return sam2bam(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2fq") == 0) return bam2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "rmdup") == 0) return bam_rmdup(argc-1, argv+1);
    else if (strcmp(argv[1], "anno") == 0) return bam_anno_attr(argc-1, argv+1);
    else if (strcmp(argv[1], "corr") == 0) return bam_corr_umi(argc-1, argv+1);
    else if (strcmp(argv[1], "attrcnt") == 0) return bam_count_attr(argc-1, argv+1);
    else if (strcmp(argv[1], "extract") == 0) return bam_extract_tags(argc-1, argv+1);
    else if (strcmp(argv[1], "pick") == 0) return bam_pick(argc-1, argv+1);
    // else if (strcmp(argv[1], "genecov") == 0) return gene_cov(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2frag") == 0) return bam2frag(argc-1, argv+1);
    else if (strcmp(argv[1], "count") == 0) return count_matrix(argc-1, argv+1);
    else if (strcmp(argv[1], "fusion") == 0) return gene_fusion(argc-1, argv+1);
    // else if (strcmp(argv[1], "assem") == 0)  return fastq_assem(argc-1, argv+1);
    // else if (strcmp(argv[1], "segment") == 0) return fastq_segment(argc-1, argv+1);
    // else if (strcmp(argv[1], "segment2") == 0) return check_segment2(argc-1, argv+1);
    //else if (strcmp(argv[1], "cleanup") == 0) return LFR_cleanup(argc-1, argv+1);
    //else if (strcmp(argv[1], "overlap") == 0) return fastq_overlap(argc-1, argv+1);
    // else if (strcmp(argv[1], "impute") == 0) return LFR_impute(argc-1, argv+1);                    
    else return usage();
    return 0;
}
