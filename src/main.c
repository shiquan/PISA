#include "utils.h"
#include "single_cell_version.h"
#include "version.h"
#include <string.h>

int usage()
{
    fprintf(stderr, "\nPISA - a Preprocessing and In silico dissecting toolkit for Single-cell data Analysis.\n");
    fprintf(stderr, "Version : %s + htslib : %s\n", SINGLECELL_VERSION, HTS_VERSION);
    fprintf(stderr, "\nCommands :\n");
    fprintf(stderr, "\n--- Processing FASTQ\n");
    fprintf(stderr, "    parse      Parse barcodes and UMIs from fastq reads.\n");
    // fprintf(stderr, "    trim       Trim TN5 mosic ends or polyAs.\n");
    fprintf(stderr, "    fsort      Sort fastq records by tags.\n");
    fprintf(stderr, "    unitig     Construct unitigs by tags.\n");
    fprintf(stderr, "    segment    Trim predefined segments from untigs.\n");
    
    fprintf(stderr, "\n--- Processing BAM\n");
    fprintf(stderr, "    sam2bam    Convert SAM to BAM, and transform tags from read name to standard SAM tags.\n");
    fprintf(stderr, "    rmdup      Remove PCR duplicates, consider cell barcodes and UMI tags.\n");
    fprintf(stderr, "    pick       Pick alignments of specific cells.\n");
    fprintf(stderr, "    anno       Annotate peak or gene names into BAM attributions.\n");
    fprintf(stderr, "    corr       Correct low frequency UMI to close one, 1 mismatch considered.\n");
    fprintf(stderr, "    attrcnt    Count raw reads and reads with predefined tag for each cell.\n");
    fprintf(stderr, "    extract    Extract tags information from BAM.\n");
    fprintf(stderr, "    count      Count matrix.\n");
    fprintf(stderr, "    bam2fq     Convert bam to fastq file. Select tags will be export at read name.\n");
    fprintf(stderr, "    genecov    Calculate reads coverage over gene body.\n");
    fprintf(stderr, "    bam2frag   Convert bam to bgzipped fragment file.\n");
    fprintf(stderr, "    impute     Impute tags in BAM.\n");
    
    // fprintf(stderr, "\n--- Processing scLFR reads. *experimental*\n");
    // fprintf(stderr, "    assem      Assem reads per fastq block (specified with -tag, fastq need be sorted)\n");
    // fprintf(stderr, "    segment    Check predefined segments from reads.\n");
    //fprintf(stderr, "    cleanup    Clean up reads before assembly.\n");   
    fprintf(stderr, "\n");
    // fprintf(stderr, "Author : Shi Quan [shiquan(AT)genomics.cn]\n");
    // fprintf(stderr, "Homepage : https://github.com/shiquan/SingleCellTools\n");
    return 1;
}
int main(int argc, char *argv[])
{
    // process FQ
    extern int fastq_prase_barcodes(int argc, char *argv[]);
    extern int fastq_trim_adaptors(int argc, char *argv[]);
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
    extern int gene_cov(int argc, char **argv);
    extern int bam2frag(int argc, char **argv);
    // process scLFR
    //extern int assem(int argc, char **argv);
    extern int check_segment(int argc, char **argv);
    extern int check_segment2(int argc, char **argv);
    extern int LFR_cleanup(int argc, char **argv);
    extern int fastq_unitig(int argc, char **argv);
    extern int LFR_impute(int argc, char **argv);

    if (argc == 1) return usage();
    else if (strcmp(argv[1], "parse") == 0) return fastq_prase_barcodes(argc-1, argv+1);
    else if (strcmp(argv[1], "trim") == 0) return fastq_trim_adaptors(argc-1, argv+1);
    else if (strcmp(argv[1], "fsort") == 0) return fsort(argc-1, argv+1);
    else if (strcmp(argv[1], "sam2bam") == 0) return sam2bam(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2fq") == 0) return bam2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "rmdup") == 0) return bam_rmdup(argc-1, argv+1);
    else if (strcmp(argv[1], "anno") == 0) return bam_anno_attr(argc-1, argv+1);
    else if (strcmp(argv[1], "corr") == 0) return bam_corr_umi(argc-1, argv+1);
    else if (strcmp(argv[1], "attrcnt") == 0) return bam_count_attr(argc-1, argv+1);
    else if (strcmp(argv[1], "extract") == 0) return bam_extract_tags(argc-1, argv+1);
    else if (strcmp(argv[1], "pick") == 0) return bam_pick(argc-1, argv+1);
    else if (strcmp(argv[1], "genecov") == 0) return gene_cov(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2frag") == 0) return bam2frag(argc-1, argv+1);
    else if (strcmp(argv[1], "count") == 0) return count_matrix(argc-1, argv+1);
    //else if (strcmp(argv[1], "assem") == 0)  return assem(argc-1, argv+1);
    else if (strcmp(argv[1], "segment") == 0) return check_segment2(argc-1, argv+1);
    // else if (strcmp(argv[1], "segment2") == 0) return check_segment2(argc-1, argv+1);
    //else if (strcmp(argv[1], "cleanup") == 0) return LFR_cleanup(argc-1, argv+1);
    else if (strcmp(argv[1], "unitig") == 0) return fastq_unitig(argc-1, argv+1);
    else if (strcmp(argv[1], "impute") == 0) return LFR_impute(argc-1, argv+1);                    
    else return usage();
    return 0;
}
