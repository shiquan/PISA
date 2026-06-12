#include <stdio.h>

int frag_count_usage()
{
    fprintf(stderr, "# Generate peak-by-cell fragment count matrix.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count2 -bed peaks.bed -t 10 -list barcodes.txt -outdir exp fragments.tsv.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -list     [FILE]     Barcode whitelist, used as column names in the matrix. If not set, all barcodes are counted.\n");
    fprintf(stderr, " -bed      [BED]      Peak regions.\n");
    fprintf(stderr, " -outdir   [DIR]      Output directory for MEX-format matrix.\n");
    fprintf(stderr, " -prefix   [STR]      Prefix for output files.\n");
    fprintf(stderr, " -t        [INT]      Threads.\n");
    fprintf(stderr, "\n");

    return 1;
}
int fragment_usage()
{
    fprintf(stderr, "# Convert SAM records to fragment file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bam2frag -cb CB -list cell_barcodes.txt -o out.tsv.gz in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o       [FILE]    Output file, bgzipped and indexed.\n");
    fprintf(stderr, " -cb      [TAG]     Cell barcode tag name.\n");
    fprintf(stderr, " -list    [FILE]    Cell barcode whitelist.\n");
    fprintf(stderr, " -q       [20]      Minimum mapping quality to filter reads.\n");
    fprintf(stderr, " -isize   [2000]    Skip fragments with insert size larger than this. [2KB]\n");
    fprintf(stderr, " -bed     [BED]     Only output fragments overlapping target regions.\n");
    fprintf(stderr, " -black-region [BED] Exclude fragments overlapping blacklist regions.\n");
    fprintf(stderr, " -stat    [FILE]    Output transposition events per cell.\n");
    fprintf(stderr, " -@       [4]       Threads for compression. [4]\n");
    fprintf(stderr, " -disable-offset    Disable Tn5 offset adjustment.\n");
    fprintf(stderr, "\n");
    return 1;
}

int fastq_parse_usage()
{
    fprintf(stderr, "# Parse cell barcode and UMI from raw FASTQ.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse -config read_struct.json -report fastq.csv -cbdis cell_dist.tsv \\\n");
    fprintf(stderr, "          -1 out.fq lane1_1.fq.gz,lane02_1.fq.gz  lane1_2.fq.gz,lane2_2.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -1       [fastq]   Read 1 output. [stdout]\n");
    fprintf(stderr, " -2       [fastq]   Read 2 output.\n");
    fprintf(stderr, " -config  [json]    Read structure configuration file in JSON format. Required.\n");
    //fprintf(stderr, " -rule    [STRING]  Read structure definition. See \x1b[31m\x1b[1mNotice\x1b[0m.\n");
    fprintf(stderr, " -run     [string]  Platform preset code for specific library kit.\n");
    fprintf(stderr, " -cbdis   [FILE]    Output read count per cell barcode.\n");
    fprintf(stderr, " -p                 Read 1 and read 2 are interleaved in the input file.\n");
    fprintf(stderr, " -q       [INT]     Drop reads with average sequencing quality below this value.\n");
    fprintf(stderr, " -dropN             Drop reads with N bases in sequence or barcode.\n");
    fprintf(stderr, " -report  [csv]     Summary report.\n");
    fprintf(stderr, " -t       [INT]     Threads. [4]\n");
    fprintf(stderr, " -dist    1|2|3     Distance scoring method: 1=Hamming, 2=Levenshtein, 3=both.\n");
    fprintf(stderr, "                    With 3, Hamming is tried first; on no match, Levenshtein is used.\n");

    //fprintf(stderr, " -x                 Preset read structure. Use one of codes predefined below.\n");
    //fprintf(stderr, "        - C4v1      MGI DNBelab C4 RNA v1/v2 kit\n");
    //fprintf(stderr, "        - 10Xv3     10X Genomics 3' v3 kit. Use barcode whitelist from \"3M-febrary-2018.txt.gz\"\n");
    //fprintf(stderr, "        - 10Xv2     10X Genomics 3' v2 or 5' v1/v2 kit. Use barcode whitelist from \"737k-august-2016.txt\"\n");
    //fprintf(stderr, "        - 10Xv3     10X Genomics 3' v1 kit. Use barcode whitelist from \"737k-april-2014_rc.txt\"\n");
    //fprintf(stderr, "        - 10XMulv1  10X Genomics Multiome (ATAC+GEX) v1 kit. Use barcode whitelist from \"737k-arc-v1.txt.gz\"\n");
    //fprintf(stderr, "        - 10XLTv1   10X Genomics 3' low throughput (LT) kit. Use barcode whitelist from \"9K-LT-march-2021.txt.gz\"\n");
    /*
    */
    return 1;
}

int fastq_parse2_usage()
{
    fprintf(stderr, "# Parse cell barcode and UMI from raw FASTQ.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse -rule CB,R1:1-10,whitelist.txt,CB,1;R1,R1:11-60;R2,R2 -report fastq.csv \\\n");
    fprintf(stderr, "           lane1_1.fq.gz,lane02_1.fq.gz  lane1_2.fq.gz,lane2_2.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -1       [fastq]   Read 1 output.\n");
    fprintf(stderr, " -2       [fastq]   Read 2 output.\n");
    fprintf(stderr, " -rule    [STRING]  Read structure definition. See \x1b[31m\x1b[1mNotice\x1b[0m.\n");
    fprintf(stderr, " -p                 Read 1 and read 2 are interleaved in the input file.\n");
    fprintf(stderr, " -q       [INT]     Drop reads with average sequencing quality below this value.\n");
    fprintf(stderr, " -dropN             Drop reads with N bases in sequence or barcode.\n");
    fprintf(stderr, " -report  [csv]     Summary report.\n");
    fprintf(stderr, " -order             Preserve input read order.\n");
    fprintf(stderr, " -t       [INT]     Threads. [4]\n");
    //fprintf(stderr, " -suffix  [STRING]  Suffix string for corrected barcode.\n");
    fprintf(stderr, " -x                 Predefined preset for a specific library kit.\n");
    fprintf(stderr, "          * C4      Library structure for DNBelab C4 RNA kit v1.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * -rule accepts a rule string in the format \"TAG,location,whitelist,corrected TAG,mismatch\".\n");
    fprintf(stderr, "   Each rule's location uses the format R[12]:start-end. Both start and end are 1-based.\n");
    fprintf(stderr, "   TAG and location are mandatory; whitelist, corrected TAG, and mismatch are optional.\n");
    fprintf(stderr, "   Multiple rules are separated by ';'. R1/R2 in location refer to raw read files;\n");
    fprintf(stderr, "   R1/R2 in the TAG refer to output read files. Examples:\n");
    fprintf(stderr, " * PISA reads FASTQ using libz, which does not support parallel I/O.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse -rule '\x1b[32mCR,R1:1-18,barcodes.txt,CB,1;\x1b[33mUR,R1:19-30;\x1b[34mR1,R2:1-100\x1b[0m' -1 read_1.fq raw_read_1.fq raw_read_2.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "#\x1b[32m CR,R1:1-18,barcodes.txt,CB,1 \x1b[0m - CR tag at read 1 positions 1-18, corrected against\n");
    fprintf(stderr, "#   barcodes.txt (Hamming distance <= 1). Corrected cell barcode stored in CB tag.\n");
    fprintf(stderr, "#\x1b[33m UR,R1:19-30\x1b[0m - UR tag at read 1 positions 19-30.\n");
    fprintf(stderr, "#\x1b[34m R1,R2:1-100\x1b[0m - Read 2 positions 1-100 are output as read 1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse -rule '\x1b[32mCR,R1:1-10,bc1.txt,CB,1;\x1b[33mCR,R1:11-20,bc2.txt,CB,1;\x1b[34mR1,R2:1-100\x1b[0m' -1 read_1.fq raw_read_1.fq raw_read_2.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "#\x1b[32m CR,R1:1-10,bc1.txt,CB,1;\x1b[33mCR,R1:11-20,bc2.txt,CB,1\x1b[0m - Cell barcode consists of two\n");
    fprintf(stderr, "#   segments: first segment at read 1 positions 1-10 (whitelist bc1.txt),\n");
    fprintf(stderr, "#   second segment at positions 11-20 (whitelist bc2.txt). Both segments are\n");
    fprintf(stderr, "#   combined after correction since they share the same corrected tag.\n");
    fprintf(stderr, "\n");
    return 1;
}
int fsort_usage()
{
    fprintf(stderr, "# Sort reads by tag values.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m fsort -tags CB,UR -list cell_barcodes_top10K.txt -@ 5 -o sorted.fq.gz in.fq\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags [TAGs]   Tags to sort by, e.g. CB,UR. Tag order matters.\n");
    //fprintf(stderr, " -dedup              Remove duplicate reads with identical tags, keeping the highest quality copy.\n");
    //fprintf(stderr, " -dup-tag [TAG]      Tag name for duplication count. Used with -dedup only. [DU]\n");
    //  fprintf(stderr, " -list    [FILE]     Whitelist for the first tag, typically cell barcodes.\n");
    fprintf(stderr, " -@       [INT]      Threads for compression.\n");
    fprintf(stderr, " -o       [fq.gz]    BGZF-compressed output FASTQ file.\n");
    fprintf(stderr, " -m       [mem]      Memory per thread. [1G]\n");
    fprintf(stderr, " -p                  Input FASTQ is interleaved paired-end.\n");
    fprintf(stderr, " -T       [prefix]   Prefix for temporary files.\n");
    // fprintf(stderr, " -report  [csv]      Summary report.\n");
    fprintf(stderr, "\n");
    return 1;
}
int fastq_stream_usage()
{
    fprintf(stderr, "# Run a user-defined script on each FASTQ+ block.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m stream -script run.sh reads.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags [TAGs]   Tags defining read blocks.\n");
    fprintf(stderr, " -script  [FILE]     User-defined shell script; processes $FQ, writes results to stdout.\n");
    fprintf(stderr, " -min     [INT]      Minimum reads per block to process. [2]\n");
    fprintf(stderr, " -keep               Output unprocessed FASTQ+ records as well.\n");
    //fprintf(stderr, " -max     [INT]      Maximum reads per block; larger blocks are downsampled. [8000]\n");
    fprintf(stderr, " -fa                 Stream FASTA instead of FASTQ.\n");
    // fprintf(stderr, " -rename            Rename output reads\n");
    fprintf(stderr, " -tmpdir  [DIR]      Temporary directory for block files.\n");
    fprintf(stderr, " -t       [INT]      Threads.\n");
    fprintf(stderr, " -o       [FILE]     Output file.\n");
    fprintf(stderr, " -nw                 Suppress warning messages.\n");
    fprintf(stderr, "\n");
    return 1;
}

/*
int assemble_usage()
{
    fprintf(stderr, "* Assemble reads from the same barcode.\n");
    fprintf(stderr, "assem [options] in.fq\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -t       [INT]      Threads.\n");
    fprintf(stderr, " -o       [fastq]    Output fastq.\n");
    fprintf(stderr, " -tag     [TAGs]     Tags of read block.\n");
    fprintf(stderr, " -p                  Input fastq is smart paired.\n");
    fprintf(stderr, " -dis     [FILE]     Assembled length distribution.\n");
    //fprintf(stderr, " -report  [csv]      Summary information.\n");
    fprintf(stderr, "\n");
    return 1;
}

int segment_usage()
{
    fprintf(stderr, "* Select defined segments from untigs.\n");
    fprintf(stderr, "Segment [options] in.fq\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, "-config   [json]      Configure file.\n");
    fprintf(stderr, "-o        [fastq]     Trimed fastq.\n");
    fprintf(stderr, "-sl       [INT]       Seed length for mapping consensus sequence.\n");
    fprintf(stderr, "-t        [INT]       Threads.\n");
    fprintf(stderr, "-tag      [TAGs]      Tags for each read block.\n");
    fprintf(stderr, "-pb       [TAG]       Phase block tag. Default tag is PB.\n");
    fprintf(stderr, "-k                    Keep all reads even no segments detected.\n");
    fprintf(stderr, "-sum      [csv]       Summary report.\n");
    fprintf(stderr, "\n");
    return 1;
}
*/
int sam2bam_usage()
{
    fprintf(stderr, "# Parse FASTQ+ read name and convert SAM to BAM.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m sam2bam -report alignment.csv -@ 5 -adjust-mapq -gtf genes.gtf -o aln.bam in.sam[.gz]\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o       [BAM]       Output file. [stdout]\n");
    fprintf(stderr, " -t       [INT]       Worker threads.\n");
    fprintf(stderr, " -mito    [string]    Mitochondrial chromosome name, for computing mitochondrial read ratio.\n");
    fprintf(stderr, " -maln    [BAM]       Write mitochondrial reads to this file (excluded from main output).\n");
    fprintf(stderr, " -@       [INT]       Threads for BAM compression.\n");
    fprintf(stderr, " -report  [csv]       Alignment summary report.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note :\n");
    fprintf(stderr, "* Reads mapped to multiple loci are typically marked as low quality and filtered\n");
    fprintf(stderr, "  downstream. For RNA-seq, if a read maps to one exonic locus and one or more\n");
    fprintf(stderr, "  non-exonic loci, the exonic locus can be prioritized as the primary alignment\n");
    fprintf(stderr, "  with mapping quality adjusted to 255. Tag MM:i:1 is added to such records.\n");
    fprintf(stderr, "* Input SAM must be sorted by read name, and the aligner should output all hits\n");
    fprintf(stderr, "  for each read.\n");
    fprintf(stderr, " -adjust-mapq         Enable mapping quality adjustment (requires -gtf).\n");
    fprintf(stderr, " -gtf     [GTF]       GTF annotation file, required for detecting exonic regions.\n");
    fprintf(stderr, " -qual    [255]       Adjusted mapping quality value.\n");
    fprintf(stderr, "\n");
    return 1;
}

int rmdup_usage()
{
    fprintf(stderr, "# Remove PCR duplicates sharing identical tag values.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m rmdup -tags CB,UR -o rmdup.bam in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags [TAGs]    Tags to group reads for deduplication.\n");
    fprintf(stderr, " -@        [INT]      Threads for BAM decompression.\n");
    fprintf(stderr, " -o        [BAM]      Output BAM file.\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold.\n");
    // fprintf(stderr, " -S                   Treat paired-end reads as single-end.\n");
    fprintf(stderr, " -k                   Keep duplicates, flag them instead of removing.\n");
    fprintf(stderr, " -nw                  Suppress warning messages.\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, "* Currently only single-end reads are supported.\n");
    fprintf(stderr, "\n");
    return 1;
}

int pick_usage()
{
    fprintf(stderr, "# Filter alignment records by barcode whitelist.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m pick -tags CB,GN -list cell_barcodes.txt in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags [TAGs]    Tag names to filter on.\n");
    fprintf(stderr, " -list     [FILE]     Barcode whitelist; values in the corresponding column are filtered.\n");
    fprintf(stderr, " -o        [BAM]      Output BAM file.\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold.\n");
    fprintf(stderr, " -@        [INT]      Threads for BAM decompression.\n");
    fprintf(stderr, "\n");
    return 1;
}
int anno_usage()
{
    fprintf(stderr, "# Annotate SAM/BAM records with overlapping functional regions.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m anno -bed peak.bed -tag PK -vcf in.vcf.gz -vtag VF -vcf-ss -ref-alt -o anno.bam in.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m anno -gtf genes.gtf -o anno.bam in.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m anno -gtf genes.gtf -o anno.bam -sam in.sam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o        [BAM]       Output BAM file.\n");
    fprintf(stderr, " -report   [csv]       Summary report.\n");
    fprintf(stderr, " -@        [INT]       Threads for BAM compression.\n");
    fprintf(stderr, " -q        [0]         Minimum mapping quality; reads below this are skipped.\n");
    fprintf(stderr, " -t        [INT]       Worker threads for annotation.\n");
    fprintf(stderr, " -chunk    [INT]       Chunk size per thread.\n");
    fprintf(stderr, " -anno-only            Export annotated reads only.\n");
    fprintf(stderr, " -sam                  Input is SAM format; parse tags from read name.\n");
    fprintf(stderr, " -rev                  Annotate on reverse strand (for probe-ligation FFPE libraries).\n");
    fprintf(stderr, " -is                   Disable strand-specific annotation.\n");
    fprintf(stderr, "\nOptions for BED file :\n");
    fprintf(stderr, " -bed      [BED]       Functional regions. 3- or 4-column BED; col 4 may be empty or contain region names.\n");
    fprintf(stderr, " -tag      [TAG]       Annotation tag name for BED regions. [PK]\n");

    fprintf(stderr, "\nOptions for mixed-species samples :\n");
    fprintf(stderr, " -chr-species  [FILE]  Chromosome-to-species mapping file.\n");
    fprintf(stderr, " -btag     [TAG]       Species tag name. [SP]\n");

    fprintf(stderr, "\nOptions for GTF file :\n");
    fprintf(stderr, " -gtf      [GTF]       GTF annotation file. gene_id and transcript_id are required for each record.\n");
    fprintf(stderr, " -tags     [TAGs]      Annotation tag names; see \x1b[31m\x1b[1mNotice\x1b[0m below. [TX,GN,GX,RE,EX,JC]\n");
    fprintf(stderr, " -splice               Annotate reads spanning exon-intron junctions with all tags.\n");
    fprintf(stderr, " -intron/-velo         Annotate reads covering intronic regions with all tags.\n");
    fprintf(stderr, " -exon                 Generate exon-level and junction annotations.\n");
    fprintf(stderr, "                       Exon name (chr:start-end/[[+-]/Gene) stored in EX tag,\n");
    fprintf(stderr, "                       junction name (chr:exon1_end-exon2_start/[+-]/Gene) in JC tag.\n");
    fprintf(stderr, " -flatten              Split overlapping exons into non-overlapping bins.\n");
    fprintf(stderr, " -psi                  Annotate excluded reads (ER tag) for each exon.\n");
    fprintf(stderr, " -vague-edge [INT]     Allow junction reads not precisely aligned to splice sites or gene boundaries.\n");
    fprintf(stderr, "                       Intended for long-read (third-generation) data.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " -as                   Annotate antisense reads; gene name saved in AT tag.\n");
    fprintf(stderr, " -tss                  Annotate reads starting at TSS, for capped RNA libraries. **experimental**\n");
    fprintf(stderr, " -ctag     [TAG]       Tag name for TSS annotation. Requires -tss.\n");

    fprintf(stderr, "\nOptions for VCF file :\n");
    fprintf(stderr, " -vcf      [VCF/BCF]   Variant file. By default, only alternative alleles are annotated.\n");
    fprintf(stderr, " -vtag     [TAG]       Tag name for variant annotation. [VR]\n");
    fprintf(stderr, " -ref-alt              Also annotate reference alleles.\n");
    // fprintf(stderr, " -vcf-ss               Annotate variants in strand sensitive way.\n");
    //fprintf(stderr, " -phased               Annotate phase block for phased positions.\n");

    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * For SAM input, tags in the read name are parsed automatically.\n");
    fprintf(stderr, " * For GTF mode, default tag names are set below; override them with -tags.\n");
    fprintf(stderr, "   TX : Transcript ID.\n");
    //fprintf(stderr, "   AN : Same with TX but set only if read mapped to antisense strand of transcript.\n");
    fprintf(stderr, "   GN : Gene name.\n");
    fprintf(stderr, "   GX : Gene ID.\n");
    fprintf(stderr, "   RE : Region type:\n");
    fprintf(stderr, "        E  = exon,  N  = intron,  C  = exon+intron overlap,\n");
    fprintf(stderr, "        S  = junction (spans isoforms),  V  = ambiguous,\n");
    fprintf(stderr, "        I  = intergenic,  A  = antisense/antisense exon,\n");
    fprintf(stderr, "        B  = antisense intron,  X  = one or more exons excluded in transcript.\n");
    fprintf(stderr, "   AT : Gene name on antisense strand (only with -as).\n");
    fprintf(stderr, " * Tags set with -exon :\n");
    fprintf(stderr, "   EX : Exon name.\n");
    fprintf(stderr, "   JC : Isoform junction name.\n");
    fprintf(stderr, "   FL : Flattened exon name (only with -flatten).\n");
    fprintf(stderr, " * Tags set with -psi :\n");
    fprintf(stderr, "   ER : Excluded exons.\n");
    fprintf(stderr, " * PSI = EX / (EX + ER), where EX counts inclusion reads, ER counts exclusion reads.\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_corr_usage()
{
    fprintf(stderr, "# Correct similar barcodes within Hamming distance 1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m corr -tag UR -new-tag UB -tags-block CB,GN -@ 5 -o corr.bam in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o        [BAM]       Output BAM file.\n");
    fprintf(stderr, " -tag      [TAG]       Tag to correct.\n");
    fprintf(stderr, " -new-tag  [TAG]       Tag name for corrected values.\n");
    fprintf(stderr, " -tags-block  [TAGs]   Tags defining read groups for correction scope.\n");
    fprintf(stderr, " -cr                   Enable Cell Ranger-style UMI correction. See \x1b[31m\x1b[1mExamples\x1b[0m below.\n");
    fprintf(stderr, " -e        [INT]       Maximum Hamming distance for matching. [1]\n");
    fprintf(stderr, " -@        [INT]       Threads for BAM compression.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mExamples\x1b[0m :\n");
    fprintf(stderr, " // Two read groups share the same cell barcode (CB) and gene (GN), but their raw\n");
    fprintf(stderr, " // UMIs (UR) differ by one base. The less-supported UMI is corrected to match the\n");
    fprintf(stderr, " // higher-support UMI. The corrected UMI is stored in tag UB.\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m corr -tag UR -new-tag UB -tags-block CB,GN in.bam -o corr.bam\n\n");
    fprintf(stderr, " // Same as above. Additionally, if multiple groups share the same CB and UB but\n");
    fprintf(stderr, " // different GN, only the GN with the most reads is kept; all others are discarded.\n");
    fprintf(stderr, " // In case of a tie, all groups are discarded (Cell Ranger method).\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m corr -cr -tag UR -new-tag UB -tags-block CB,GN in.bam -o corr.bam\n\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_attr_usage()
{
    fprintf(stderr, "# Count the frequency of tag values.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m attrcnt -cb CB -tag UR,GN -dedup -all-tags in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -cb       [TAG]      Cell barcode or other tag for grouping reads.\n");
    fprintf(stderr, " -list     [FILE]     Cell barcode whitelist.\n");
    fprintf(stderr, " -tag/-tags [TAGs]    Tags to count.\n");
    fprintf(stderr, " -dedup               Deduplicate attributes within each tag.\n");
    fprintf(stderr, " -all-tags            Only count records that have all specified tags.\n");
    fprintf(stderr, " -group    [TAG]      Group tag; count tags separately for each group.\n");
    fprintf(stderr, " -o        [FILE]     Output count table.\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold.\n");
    fprintf(stderr, " -no-header           Omit header from output.\n");
    fprintf(stderr, " -@        [INT]      Threads for BAM decompression.\n");
    fprintf(stderr, " -ttag     [TAG]      Region type tag. [RE]\n");
    fprintf(stderr, " -ttype    [TYPE]     Region types to count. E.g. 'E,S' for exonic reads, 'N,C' for intronic.\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_extract_usage()
{
    fprintf(stderr, "# Extract tag values from alignments.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m extract -tag CB,UR,GN -o tags.tsv in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags [TAGs]    Tags to extract.\n");
    fprintf(stderr, " -o        [FILE]     Output file in TSV format.\n");
    fprintf(stderr, " -n                   Include read name in output.\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold.\n");
    fprintf(stderr, " -all                 Only export records where all specified tags have values.\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_count_usage()
{
    fprintf(stderr, "# Generate count matrix for single-cell datasets.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count -cb CB -anno-tag GN -umi UB -outdir exp aln.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count [options] aln1.bam,aln2.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count -cb RG -sample-list bam_files.txt -outdir exp\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count -tags Cx,Cy -anno-tag GN -umi UB -outdir exp -velo aln.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags/-cb [TAGs] Cell barcode tag, or two coordinate tags for spatial data.\n");
    fprintf(stderr, " -anno-tag [TAGs]     Annotation tag (gene or peak). With multiple tags, the first\n");
    fprintf(stderr, "                      is counted only when all other tags are present.\n");
    fprintf(stderr, " -genome-bin [INT]    Generate genome bin count matrix. Conflicts with -anno-tag and -chr.\n");
    fprintf(stderr, " -is                  Ignore strand for bin counting.\n");
    fprintf(stderr, " -chr                 Generate per-chromosome count. Conflicts with -anno-tag and -genome-bin.\n");
    fprintf(stderr, " -list     [FILE]     Barcode whitelist. If not set, all barcodes are counted.\n");
    fprintf(stderr, " -outdir   [DIR]      Output directory for MEX-format matrix.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. Each UMI is counted once per gene or peak.\n");
    fprintf(stderr, " -one-hit             Skip reads that map to more than one gene or peak.\n");
    // fprintf(stderr, " -corr                Enable UMI correction (Hamming distance <= 1).\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold. [20]\n");
    fprintf(stderr, " -t        [INT]      Threads.\n");
    fprintf(stderr, " -ttag     [TAG]      Region type tag. [RE]\n");
    fprintf(stderr, " -velo                Generate spliced and unspliced matrices for RNA velocity.\n");
    fprintf(stderr, " -ttype    [TYPE]     Region types to count. E.g. 'E,S' for exonic, 'N,C' for intronic.\n");
    //fprintf(stderr, " -file-barcode        Use BAM file name as cell barcode. Requires -sample-list.\n");
    fprintf(stderr, " -sample-list [FILE]  File listing BAM paths, one per line.\n");

    fprintf(stderr, "\nOptions for Stereo-seq :\n");
    fprintf(stderr, " -stereoseq           Stereo-seq UMI encoding (hex string); enables decoding.\n");
    fprintf(stderr, " -spatial-bin [INT]   Bin size for spatial coordinates. [1]\n");
    fprintf(stderr, " -dup                 Do not skip duplicate reads.\n");

    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * Region type (RE), annotated by `\x1b[1mPISA\x1b[0m anno`, labels the functional region.\n");
    fprintf(stderr, "   Use -ttype with one or more types (E/S/C/N) to count only those regions.\n");
    fprintf(stderr, " * To count from multiple BAM files, separate paths with ',' or use -sample-list.\n");
    // fprintf(stderr, " * -cb conflicts with -file-barcode.\n");
    fprintf(stderr, " * With -velo, spliced/ and unspliced/ folders are created in -outdir.\n");
    fprintf(stderr,"\n");
    return 1;
}

int bam2fq_usage()
{
    fprintf(stderr, "# Convert BAM to FASTQ+.\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bam2fq -tag CB,UB,GN -o out.fq aln.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag/-tags [TAGs]    Tags to embed in output read names.\n");
    fprintf(stderr, " -f                   Skip records missing any of the specified tags.\n");
    fprintf(stderr, " -fa                  Output FASTA instead of FASTQ.\n");
    fprintf(stderr, " -o        [fastq]    Output file.\n");
    fprintf(stderr, " -@        [INT]      Threads for BAM decompression.\n");
    fprintf(stderr, "\n");
    //fprintf(stderr, "* Following options are experimental.\n");
    //fprintf(stderr, "* Merge overlapping reads from the same molecule.\n");
    //fprintf(stderr, " -tag      [TAGs]     Tags to group reads.\n");
    //fprintf(stderr, " -i        [TAGs]     Tags to inherit on merged reads.\n");
    //fprintf(stderr, " -strand              Only merge reads on the same strand.\n");
    //fprintf(stderr, "\n");
    return 1;
}
/*
int bam_impute_usage()
{
    fprintf(stderr, "* Imputate empty tag by existed tags.\n");
    fprintf(stderr, "bam_impute in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, "  -impute  [TAGs]     Tags to impute.\n");
    fprintf(stderr, "  -block   [TAGs]     Tags to identify each block.\n");
    fprintf(stderr, "  -dist    [INT]      Distance between reads from same block will be imputed.\n");
    fprintf(stderr, "  -k                   Keep unclassified reads in the output.\n");
    fprintf(stderr, "  -@       [INT]      Threads to unpack BAM.\n");
    fprintf(stderr, "\n");
    return 1;
}
*/

int gene_fusion_usage()
{
    fprintf(stderr, "# Predict gene fusion events. **experimental**\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m fusion [options] in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -gn    [TAG]     Gene name tag. [GN]\n");
    fprintf(stderr, " -cb    [TAG]     Cell barcode tag. [CB]\n");
    fprintf(stderr, " -umi   [TAG]     UMI tag. [UB]\n");
    fprintf(stderr, " -fs    [TAG]     Fusion tag name written to output. [FS]\n");
    fprintf(stderr, " -o     [BAM]     Output BAM file.\n");
    fprintf(stderr, " -q     [INT]     Minimum mapping quality threshold. [0]\n");
    fprintf(stderr, " -list  [FILE]    Cell barcode whitelist.\n");
    fprintf(stderr, " -@     [INT]     Threads for BAM decompression.\n");
    fprintf(stderr, " -m     [INT]     Maximum genes per UMI. [3]\n");
    fprintf(stderr, " -fusion-only      Export only fusion-supporting reads.\n");
    fprintf(stderr, "\n");
    return 1;
}

int depth_usage()
{
    fprintf(stderr, "# Compute coverage depth or unique UMI count per genomic position.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m depth sorted.bam chr1:1-2\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m depth sorted.bam chr1:1-2:+\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m depth -tag CB -list cells.txt -umi UB -region in.bed -o depth.tsv sorted.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m depth -tag CB -list cells.txt -umi UB sorted.bam chr1:1-2:+\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag      [TAG]      Tag for grouping reads.\n");
    fprintf(stderr, " -list     [FILE]     Whitelist for -tag.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. When set, only unique UMIs are counted per position.\n");
    fprintf(stderr, " -split               Group reads by tag value.\n");
    fprintf(stderr, " -is                  Ignore strand.\n");
    fprintf(stderr, " -bed      [BED]      Target BED regions.\n");
    fprintf(stderr, " -o        [FILE]     Output file. [stdout]\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold. [20]\n");
    fprintf(stderr, " -@        [INT]      Threads for BAM decompression. [4]\n");
    fprintf(stderr, " -0                   Also output positions with zero coverage.\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * Requires a coordinate-sorted and indexed BAM.\n");
    fprintf(stderr, " * Unlike samtools depth, PISA depth accounts for UMIs and strand.\n");
    return 1;
}

int bin_usage()
{
    fprintf(stderr, "# Count reads or fragments per genomic bin.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bin -outdir bins sorted.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bin -cb CB -outdir bins sorted.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bin -cb CB -umi UB -outdir bins sorted.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -cb       [TAG]      Cell barcode tag, or two coordinate tags for spatial data.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. Each UMI is counted once per bin.\n");
    fprintf(stderr, " -bins     [SIZES]    Comma-separated bin sizes. [\"50K,100K,500K,1M\"]\n");
    fprintf(stderr, " -list     [FILE]     Barcode whitelist. If not set, all barcodes are counted.\n");
    fprintf(stderr, " -outdir   [DIR]      Output directory for MEX-format count files.\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold. [20]\n");
    fprintf(stderr, " -strand              Enable strand-specific counting.\n");
    fprintf(stderr, " -t        [INT]      Threads.\n");
    fprintf(stderr, " -sample-list [FILE]  File listing BAM paths with sample names, one per line.\n");
    fprintf(stderr,"\n");
    return 1;
}

int cov_usage()
{
    fprintf(stderr, "# Count genomic bases with depth >= min-depth per tag.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m cov -tag CB -o cov.csv sorted.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m cov -tag CB -list cells.txt -umi UB -min-depth 5 -o cov.tsv sorted.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag      [TAG]      Barcode/tag to group reads. Required.\n");
    fprintf(stderr, " -list     [FILE]     Whitelist for -tag.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. When set, only unique UMIs count toward depth.\n");
    fprintf(stderr, " -min-depth  [INT]    Minimum depth to count a base as covered. [%d]\n", 1);
    fprintf(stderr, " -o        [FILE]     Output file (CSV or TSV, depending on extension). [stdout]\n");
    fprintf(stderr, " -q        [INT]      Minimum mapping quality threshold. [%d]\n", 20);
    fprintf(stderr, " -@        [INT]      Threads for BAM decompression. [%d]\n", 4);
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * Requires a coordinate-sorted and indexed BAM.\n");
    fprintf(stderr, " * With -umi, only unique UMIs are counted toward depth per position per cell.\n");
    fprintf(stderr, " * With -min-depth, only bases with depth >= N are counted as covered.\n");
    fprintf(stderr, " * Output format is CSV (.csv) or TSV (.tsv) depending on output filename.\n");
    fprintf(stderr, " * With -list, only reads with barcodes in the whitelist are counted.\n");
    return 1;
}
