#include <stdio.h>

int frag_count_usage()
{
    fprintf(stderr, "# Count fragments per peak per cell matrix.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count2 -bed peaks.bed -t 10 -list barcodes.txt -outdir exp fragments.tsv.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -list     [FILE]     Barcode white list, used as column names at matrix. If not set, all barcodes will be count.\n");
    fprintf(stderr, " -bed      [BED]      Peaks.\n");
    fprintf(stderr, " -outdir   [DIR]      Output matrix in MEX format into this fold.\n");
    fprintf(stderr, " -prefix   [STR]      Prefix of output files.\n");
    fprintf(stderr, " -t        [INT]      Threads.\n");
    fprintf(stderr, "\n");

    return 1;
}
int fragment_usage()
{
    fprintf(stderr, "# Convert sam record to fragment file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bam2frag -cb CB -list cell_barcodes.txt -o out.tsv.gz in.bam\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, " -o       [FILE]    Output file. This file will be bgzipped and indexed.\n");
    fprintf(stderr, " -cb      [TAG]     Cell barcode tag.\n");
    fprintf(stderr, " -list    [FILE]    Cell barcode white list.\n");
    fprintf(stderr, " -q       [20]      Mapping quality score to filter reads.\n");
    fprintf(stderr, " -isize   [2000]    Skip if insert size greater than this. [2KB]\n");
    fprintf(stderr, " -bed     [BED]     Only convert fragments overlapped with target regions.\n");
    fprintf(stderr, " -black-region [BED] Skip convert fragments overlapped with black regions.\n");
    fprintf(stderr, " -stat    [FILE]    Transposition events per cell.\n");   
    fprintf(stderr, " -@       [4]       Thread to unpack and pack files.[4]\n");
    fprintf(stderr, " -disable-offset    Disable Tn5 offset for each fragment.\n");
    fprintf(stderr, "\n");
    return 1;
}

int fastq_parse_usage()
{
    fprintf(stderr, "# Parse cell barcode and UMI string from raw FASTQ.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse -config read_struct.json -report fastq.csv -cbdis cell_dist.tsv \\\n");
    fprintf(stderr, "          -1 out.fq lane1_1.fq.gz,lane02_1.fq.gz  lane1_2.fq.gz,lane2_2.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -1       [fastq]   Read 1 output. Default is stdout.\n");
    fprintf(stderr, " -2       [fastq]   Read 2 output.\n");
    fprintf(stderr, " -config  [json]    Read structure configure file in JSON format. Required.\n");
    //fprintf(stderr, " -rule    [STRING]  Read structure in line. See \x1b[31m\x1b[1mNotice\x1b[0m.\n");
    fprintf(stderr, " -run     [string]  Run code, used for different library.\n");
    fprintf(stderr, " -cbdis   [FILE]    Read count per cell barcode.\n");
    fprintf(stderr, " -p                 Read 1 and read 2 interleaved in the input file.\n");
    fprintf(stderr, " -q       [INT]     Drop reads if average sequencing quality below this value.\n");
    fprintf(stderr, " -dropN             Drop reads if N base in sequence or barcode.\n");
    fprintf(stderr, " -report  [csv]     Summary report.\n");
    fprintf(stderr, " -t       [INT]     Threads. [4]\n");
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
    fprintf(stderr, "# Parse cell barcode and UMI string from raw FASTQ.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse2 -rule CB,R1:1-10,whitelist.txt,CB,1;R1,R1:11-60;R2,R2 -report fastq.csv \\\n");
    fprintf(stderr, "           lane1_1.fq.gz,lane02_1.fq.gz  lane1_2.fq.gz,lane2_2.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -1       [fastq]   Read 1 output.\n");
    fprintf(stderr, " -2       [fastq]   Read 2 output.\n");
    fprintf(stderr, " -rule    [STRING]  Read structure in line. See \x1b[31m\x1b[1mNotice\x1b[0m.\n");
    fprintf(stderr, " -p                 Read 1 and read 2 interleaved in the input file.\n");
    fprintf(stderr, " -q       [INT]     Drop reads if average sequencing quality below this value.\n");
    fprintf(stderr, " -dropN             Drop reads if N base in sequence or barcode.\n");
    fprintf(stderr, " -report  [csv]     Summary report.\n");
    fprintf(stderr, " -order             Keep input order.\n");
    fprintf(stderr, " -t       [INT]     Threads. [4]\n");
    fprintf(stderr, " -x                 Predefined code for specific library.\n");
    fprintf(stderr, "          * C4      Library structure for DNBelab C4 RNA kit v1.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    //fprintf(stderr, " * \x1b[1mPISA\x1b[0m parse requires -config, -rule, or -x option to specify cell barcode and UMI locations in the raw read files.\n");
    //fprintf(stderr, " * -config accept configure file in JSON format. Full details to generate this file can be found at wiki page.\n");
    fprintf(stderr, " * -rule accept tag rule STRING to parse input fastq following format \"TAG,location,whitelist,corrected TAG,allow mismatch\".\n");
    fprintf(stderr, "   For each tag rule, location part should be format like R[12]:start-end. Both start and end location start from 1.\n");
    fprintf(stderr, "   TAG and locaion parts are mandatory, and whitelist, corrected TAG and mismatch are optional.\n");
    fprintf(stderr, "   Futhermore, multiply tags seperated by \';\'. In location part, R1 stands for raw read 1, R2 stands for raw read 2.\n");
    fprintf(stderr, "   In tag part, R1 stands for output read 1 while R2 stands for output read 2. Here are some examples.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse2 -rule '\x1b[32mCR,R1:1-18,barcodes.txt,CB,1;\x1b[33mUR,R1:19-30;\x1b[34mR1,R2:1-100\x1b[0m' -1 read_1.fq raw_read_1.fq raw_read_2.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "#\x1b[32m CR,R1:1-18,barcodes.txt,CB,1 \x1b[0m - CR tag start from 1 to 18 in read 1, and barcodes.txt are barcode whitelist,\n");
    fprintf(stderr, "#   each barcode per line. Cell barcode will be corrected while hamming distance <= 1.\n");
    fprintf(stderr, "#   Corrected cell barcode tag is CB. \n");
    fprintf(stderr, "#\x1b[33m UR,R1:19-30\x1b[0m - UR tag start from 19-30 in read 1.\n");
    fprintf(stderr, "#\x1b[34m R1,R2:1-100\x1b[0m - Sequence from 1 to 100 in read 2 output to read 1 file. \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m parse2 -rule '\x1b[32mCR,R1:1-10,bc1.txt,CB,1;\x1b[33mCR,R1:11-20,bc2.txt,CB,1;\x1b[34mR1,R2:1-100\x1b[0m' -1 read_1.fq raw_read_1.fq raw_read_2.fq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "#\x1b[32m CR,R1:1-10,bc1.txt,CB,1;\x1b[33mCR,R1:11-20,bc2.txt,CB,1\x1b[0m - This cell barcode consist of two segments, first segment start\n");
    fprintf(stderr, "#   from 1 to 10 in read 1, and whitelist is bc1.txt, and second segment start from 11 to 20, and whitelist is bc2.txt.\n");
    fprintf(stderr, "#   These two segments will be combined after correction, because the corrected tag are the same.\n");
    fprintf(stderr, "\n");
    return 1;
}
int fsort_usage()
{
    fprintf(stderr, "# Sort reads by tags.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m fsort -tags CB,UR -list cell_barcodes_top10K.txt -@ 5 -o sorted.fq.gz in.fq\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, " -tags    [TAGS]     Tags, such as CB,UR. Order of these tags is sensitive.\n");
    //fprintf(stderr, " -dedup              Remove dna copies with same tags. Only keep reads have the best quality.\n");
    //fprintf(stderr, " -dup-tag [TAG]      Tag name of duplication counts. Use with -dedup only. [DU]\n");
    //  fprintf(stderr, " -list    [FILE]     White list for first tag, usually for cell barcodes.\n");
    fprintf(stderr, " -@       [INT]      Threads to compress file.\n");
    fprintf(stderr, " -o       [fq.gz]    bgzipped output fastq file.\n");
    fprintf(stderr, " -m       [mem]      Memory per thread. [1G]\n");
    fprintf(stderr, " -p                  Input fastq is smart pairing.\n");
    fprintf(stderr, " -T       [prefix]   Write temporary files to PREFIX.nnnn.tmp\n");
//     fprintf(stderr, " -report  [csv]      Summapry report.\n");
    fprintf(stderr, "\n");
    return 1;
}
int fastq_stream_usage()
{
    fprintf(stderr, "# Perform user-defined script for each FASTQ+ block.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m stream -script run.sh reads.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tags    [TAGS]     Tags to define read blocks.\n");
    fprintf(stderr, " -script  [FILE]     User defined bash script, process $FQ and generate results to stdout.\n");
    fprintf(stderr, " -min     [INT]      Mininal reads per block to process.  [2]\n");
    fprintf(stderr, " -keep               Output unprocessed FASTQ+ records.\n");
    //fprintf(stderr, " -max     [INT]      Maximal reads per block, if more reads, will downsampling. [8000]\n");
    fprintf(stderr, " -fa                 Stream FASTQ output instead of FASTQ.\n");
    // fprintf(stderr, " -rename           Rename output reads\n");
    fprintf(stderr, " -tmpdir\n");
    fprintf(stderr, " -t       [INT]      Threads.\n");
    fprintf(stderr, " -o       [FILE]     Path to output file.\n");
    fprintf(stderr, " -nw                 Disable warning messages.\n");
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
    fprintf(stderr, " -tag     [TAGS]     Tags of read block.\n");
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
    fprintf(stderr, "-tag      [TAGS]      Tags for each read block.\n");
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
    fprintf(stderr, " -o       [BAM]       Output file [stdout].\n");
    fprintf(stderr, " -t       [INT]       Work threads.\n");
    fprintf(stderr, " -mito    [string]    Mitochondria name. Used to stat ratio of mitochondria reads.\n");
    fprintf(stderr, " -maln    [BAM]       Export mitochondria reads into this file instead of standard output file.\n");
    fprintf(stderr, " -@       [INT]       Threads to compress bam file.\n");
    fprintf(stderr, " -report  [csv]       Alignment report.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note :\n");
    fprintf(stderr, "* Reads map to multiple loci usually be marked as low quality and filtered at downstream analysis.\n");
    fprintf(stderr, "  But for RNAseq library, if reads map to an exonic locus but also align to 1 or more non-exonic loci,\n");
    fprintf(stderr, "  the exonic locus can be prioritized as primary alignments, and mapping quality adjusts to 255. Tag\n");
    fprintf(stderr, "  MM:i:1 will also be added for this record. Following options used to adjust mapping quality.\n");
    fprintf(stderr, "* Input SAM need be sorted by read name, and aligner should output all hits of a read in this SAM.\n");
    fprintf(stderr, " -adjust-mapq         Enable adjusts mapping quality score.\n");
    fprintf(stderr, " -gtf     [GTF]       GTF annotation file. This file is required to check the exonic regions.\n");
    fprintf(stderr, " -qual    [255]       Updated quality score.\n");
    fprintf(stderr, "\n");
    return 1;    
}

int rmdup_usage()
{
    fprintf(stderr, "# Deduplicate PCR reads with same barcodes.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m rmdup -tags CB,UR -o rmdup.bam in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, "   -tags  [TAGS]       Barcode tags to group reads.\n");
    fprintf(stderr, "   -@     [INT]        Threads to unpack BAM.\n");
    fprintf(stderr, "   -o     [BAM]        Output bam.\n");
    fprintf(stderr, "   -q     [INT]        Map Quality Score cutoff.\n");
    // fprintf(stderr, "   -S                  Treat PE reads as SE.\n");
    fprintf(stderr, "   -k                  Keep duplicates, make flag instead of remove them.\n");
    fprintf(stderr, "   -nw                 Disable warnings.\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m:\n");
    fprintf(stderr, "* Currently only support single-end reads.\n");
    fprintf(stderr, "\n");
    return 1;
}

int pick_usage()
{
    fprintf(stderr, "# Pick alignment records with barcode list.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m pick -tags CB,GN -list cell_barcodes.txt in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tags    [TAGS]       Barcode tags.\n");
    fprintf(stderr, " -list    [FILE]       Barcode white list, tag values in related column will be apply.\n");
    fprintf(stderr, " -o       [BAM]        Output file.\n");
    fprintf(stderr, " -q       [INT]        Map Quality Score cutoff.\n");
    fprintf(stderr, " -@       [INT]        Threads to unpack BAM.\n");
    fprintf(stderr, "\n");
    return 1;
}
int anno_usage()
{
    fprintf(stderr, "# Annotate SAM/BAM records with overlapped function regions. Such as gene, transcript etc.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m anno -bed peak.bed -tag PK -vcf in.vcf.gz -vtag VF -o anno.bam in.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m anno -gtf genes.gtf -o anno.bam in.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m anno -gtf genes.gtf -o anno.bam -sam in.sam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o        [BAM]       Output bam file.\n");
    fprintf(stderr, " -report   [csv]       Summary report.\n");
    fprintf(stderr, " -@        [INT]       Threads to compress bam file.\n");
    fprintf(stderr, " -q        [0]         Map Quality Score cutoff. MapQ smaller and equal to this value will not be annotated.\n");
    fprintf(stderr, " -t        [INT]       Threads to annotate.\n");
    fprintf(stderr, " -chunk    [INT]       Chunk size per thread.\n");
    fprintf(stderr, " -anno-only            Export annotated reads only.\n");
    fprintf(stderr, " -sam                  Input is SAM file, parse tags from read name.\n");
    
    fprintf(stderr, "\nOptions for BED file :\n");
    fprintf(stderr, " -bed      [BED]       Function regions. Three or four columns bed file. Col 4 could be empty or names of this region.\n");
    fprintf(stderr, " -tag      [TAG]       Attribute tag name. Set with -bed.\n");

    fprintf(stderr, "\nOptions for mixed samples.\n");
    fprintf(stderr, " -chr-species  [FILE]  Chromosome name and related species binding list.\n");
    fprintf(stderr, " -btag     [TAG]       Species tag name. Set with -chr-species.\n");

    fprintf(stderr, "\nOptions for GTF file :\n");
    fprintf(stderr, " -gtf      [GTF]       GTF annotation file. gene_id,transcript_id is required for each record.\n");
    fprintf(stderr, " -tags     [TAGS]      Attribute names, more details see `\x1b[31m\x1b[1mNotice\x1b[0m` below. [TX,GN,GX,RE]\n");
    fprintf(stderr, " -ignore-strand        Ignore strand of transcript in GTF. Reads mapped to antisense transcripts will also be annotated.\n");
    fprintf(stderr, " -splice               Reads covered exon-intron edge will also be annotated with all tags.\n");
    fprintf(stderr, " -intron               Reads covered intron regions will also be annotated with all tags.\n");
    fprintf(stderr, " -tss                  Annotate reads start from TSS, designed for capped library. **experiment**\n");
    fprintf(stderr, " -ctag     [TAG]       Tag name for TSS annotation. Need set with -tss.\n");

    fprintf(stderr, "\nOptions for VCF file :\n");
    fprintf(stderr, " -vcf      [VCF/BCF]   Varaints.\n");
    fprintf(stderr, " -vtag     [TAG]       Tag name. Set with -vcf.\n");
    
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * If input is SAM format, will try to parse the tags in the read name.\n");
    fprintf(stderr, " * For GTF mode, this program will set tags in default, you could also reset them by -tags.\n");
    fprintf(stderr, "   TX : Transcript id.\n");
    //fprintf(stderr, "   AN : Same with TX but set only if read mapped to antisense strand of transcript.\n");
    fprintf(stderr, "   GN : Gene name.\n");
    fprintf(stderr, "   GX : Gene ID.\n");
    fprintf(stderr, "   RE : Region type, E (Exon), N (Intron), C (Exon and Intron), S (junction reads cover isoforms properly), V (ambiguous reads), I (Intergenic), A (Anitisense)\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_corr_usage()
{
    fprintf(stderr, "# Correct similar barcodes (hamming distance == 1).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m corr -tag UR -new-tag UB -tags-block CB,GN -@ 5 -o corr.bam in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o        [BAM]       Output bam.\n");
    fprintf(stderr, " -tag      [TAG]       Tag to correct.\n");
    fprintf(stderr, " -new-tag  [TAG]       Create a new tag for corrected barcodes.\n");
    fprintf(stderr, " -tags-block  [TAGS]   Tags to define read group. For example, if set to GN (gene), reads in the same gene will be grouped together.\n");
    fprintf(stderr, " -cr                   Enable CellRanger like UMI correction method. See `\x1b[31m\x1b[1mExamples\x1b[0m` for details.\n");
    fprintf(stderr, " -e                    Maximal hamming distance to define similar barcode, default is 1.\n");
    fprintf(stderr, " -@        [INT]       Thread to compress BAM file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mExamples\x1b[0m :\n");
    fprintf(stderr, " // Two groups of reads have same cell barcode (CB) and gene (GN) but their raw UMIs (UR) differ by only one base. The UMI of less \n");
    fprintf(stderr, " // supported is corrected to the UMI with higher support. UB save the checked or corrected UMI.\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m corr -tag UR -new-tag UB -tags-block CB,GN in.bam -o corr.bam \n\n");
    fprintf(stderr, " // Same with above. Besides, if two or more groups of reads have same CB and UB but different GN, the GN with the most supporting reads\n");
    fprintf(stderr, " // is kept for UMI counting, and the other read groups are discarded. In case of a tie for maximal read support, all read groups are\n");
    fprintf(stderr, " // discarded, as the gene cannot be confidently assigned (Cell Ranger method).\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m corr -cr -tag UR -new-tag UB -tags-block CB,GN in.bam -o corr.bam \n\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_attr_usage()
{
    fprintf(stderr, "# Count the frequency of tag values.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m attrcnt -cb CB -tags UR,GN -dedup -all-tags in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -cb       [TAG]      Cell Barcode, or other tag used for grouping reads.\n");
    fprintf(stderr, " -list     [FILE]     Cell barcode white list.\n");
    fprintf(stderr, " -tags     [TAGS]     Tags to count.\n");
    fprintf(stderr, " -dedup               Deduplicate the atrributes in each tag.\n");
    fprintf(stderr, " -all-tags            Only records with all tags be count.\n");
    fprintf(stderr, " -group    [TAG]      Group tag, count all tags for each group seperately.\n");
    fprintf(stderr, " -o        [FILE]     Output count table.\n");
    fprintf(stderr, " -q        [INT]      Map Quality to filter bam.\n");
    fprintf(stderr, " -no-header           Ignore header in the output.\n");
    fprintf(stderr, " -@        [INT]      Thread to unpack bam.\n");
    fprintf(stderr, " -ttag     [TAG]      Region type tag. [RE]\n");
    fprintf(stderr, " -ttype               Region type used to count. Set `E,S` to count exon enclosed reads. Set `N,C` to count intron overlapped reads.\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_extract_usage()
{
    fprintf(stderr, "# Extract tag values from alignments.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m extract -tags CB,UR,GN -o tags.tsv in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tags     [TAGS]     Tags to be extracted.\n");
    fprintf(stderr, " -o        [FILE]     Output file. tsv format\n");
    fprintf(stderr, " -n                   Print read name.\n");
    fprintf(stderr, " -q                   Map Quality Score threshold.\n");
    fprintf(stderr, " -all                 Only export if all tags have value.\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_count_usage()
{
    fprintf(stderr, "# Count reads or fragments matrix for single-cell datasets.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count -cb CB -anno-tag GN -umi UB -outdir exp aln.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count [options] aln1.bam,aln2.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m count -file-barcode -sample-list bam_files.txt -outdir exp\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -cb       [TAG]      Cell barcode tag.\n");
    fprintf(stderr, " -anno-tag [TAG]      Annotation tag, gene or peak.\n");
    fprintf(stderr, " -list     [FILE]     Barcode white list, used as column names at matrix. If not set, all barcodes will be count.\n");
    //fprintf(stderr, " -o        [FILE]     Output matrix.\n");
    fprintf(stderr, " -outdir   [DIR]      Output matrix in MEX format into this fold.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. Count once if more than one record has same UMI in one gene or peak.\n");
    fprintf(stderr, " -one-hit             Skip if a read hits more than 1 gene or peak.\n");
    // fprintf(stderr, " -corr                Enable correct UMIs. Similar UMIs defined as amming distance <= 1.\n");
    fprintf(stderr, " -q        [INT]      Minimal map quality to filter. Default is 20.\n");
    fprintf(stderr, " -@        [INT]      Threads to unpack BAM.\n");
    fprintf(stderr, " -ttag     [TAG]      Region type tag. [RE]\n");
    fprintf(stderr, " -velo                Generate spliced and unspliced matrix files for RNA velocity analysis.\n");
    fprintf(stderr, " -ttype    [TYPE]     Region type used to count. Set `E,S` to count exon enclosed reads. Set `N,C` to count intron overlapped reads.\n");
    fprintf(stderr, " -file-barcode        No cell barcode tag in the bam, but alias file name as cell barcode. This option must use with -sample-list.\n");
    fprintf(stderr, " -sample-list [FILE]  A list of bam files. First column of this file should be path of bam files. Optional second column is the \n");
    fprintf(stderr, "                      sample or cell name. This option is useful for one cell per bam experiment, like Smartseq.\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * Region type (RE), which label functional region reads mapped, is annotated by `\x1b[1mPISA\x1b[0m anno`. Optional -ttype can be set\n");
    fprintf(stderr, "   to one of region types(E/S/C/N) or combination to count reads mapped to these functional regions only.\n");
    fprintf(stderr, " * If you want count from more than one bam file, there are two ways to set the parameter. By seperating bam files with ',' or by\n");
    fprintf(stderr, "   setting -sample-list option. But if you want alias each bam with a predefined cell name, only -sample-list supported.\n");
    fprintf(stderr, " * -cb conflict with -file-barcode. \x1b[1mPISA\x1b[0m read cell barcode from bam tag or alias name list. Not both.\n");
    fprintf(stderr, " * If -velo set, spliced and unspliced folders will be created at outdir.\n");
    fprintf(stderr,"\n");
    return 1;
}

int bam2fq_usage()
{
    fprintf(stderr, "# Convert BAM into fastq.\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m bam2fq -tags CB,UB,GN -o out.fq aln.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tags     [TAGS]     Export tags in read name.\n");        
    fprintf(stderr, " -f                   Filter records if specified tags not all exist.\n");
    fprintf(stderr, " -fa                  Output fasta instead of fastq.\n");
    fprintf(stderr, " -o        [fastq]    Output file.\n");
    fprintf(stderr, " -@        [INT]      Threads to unpack BAM.\n");
    fprintf(stderr, "\n");
    //fprintf(stderr, "* Following options are experimental. \n");
    //fprintf(stderr, "* Merge overlapped reads from same molecular.\n");
    //fprintf(stderr, " -tag      [TAGS]     Tags to group reads.\n");
    //fprintf(stderr, " -i        [TAGS]     Inhert these tags to merged reads.\n");
    //fprintf(stderr, " -strand              Only check the overlapped reads from same strand.\n");
    //fprintf(stderr, "\n");
    return 1;
}
/*
int bam_impute_usage()
{
    fprintf(stderr, "* Imputate empty tag by existed tags.\n");
    fprintf(stderr, "bam_impute in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, "  -impute  [TAGS]     Tags to impute.\n");
    fprintf(stderr, "  -block   [TAGS]     Tags to identify each block.\n");
    fprintf(stderr, "  -dist    [INT]      Distance between reads from same block will be imputed.\n");
    fprintf(stderr, "  -k                   Keep unclassified reads in the output.\n");
    fprintf(stderr, "  -@       [INT]      Threads to unpack BAM.\n");
    fprintf(stderr, "\n");
    return 1;
}
*/

int gene_fusion_usage()
{
    fprintf(stderr, "# gene_fusion ** experiment **\n");
    fprintf(stderr, " -gn    Gene tag name in the bam. [GN]\n");
    fprintf(stderr, " -cb    Cell barcode tag name in the bam. [CB]\n");
    fprintf(stderr, " -umi   UMI barcode tag name in the bam. [UB]\n");
    fprintf(stderr, " -fs    Fusion record tag name for write into the bam. [FS]\n");
    fprintf(stderr, " -o     Output bam file.\n");
    fprintf(stderr, " -q     Mapping quality threshold. [0]\n");
    fprintf(stderr, " -list  Cell barcode white list.\n");
    fprintf(stderr, " -@     Threads to unpack bam.\n");
    fprintf(stderr, " -m     Maximal mapped genes per UMI. [3]\n");
    fprintf(stderr, " -fusion-only  Export fusion support reads only.\n");
    fprintf(stderr, "\n");
    return 1;
}

int depth_usage()
{
    fprintf(stderr, "# Count coverage depth or unique UMIs for genome locations.\n");
    fprintf(stderr, "\nUsage : \x1b[1mPISA\x1b[0m depth [options] sorted.bam [region]\n\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m depth -cb CB -umi UB -tags GN -region in.bed -o depth.tsv sorted.bam\n");
    fprintf(stderr, "\x1b[36m\x1b[1m$\x1b[0m \x1b[1mPISA\x1b[0m depth -cb CB -umi UB sorted.bam chr1:1-2:+\n");
    fprintf(stderr, "\nOptions : \n");
    fprintf(stderr, " -tag      [TAG]      Tag used for grouping reads.\n");
    fprintf(stderr, " -list     [FILE]     Candidate list for -tag.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. If set, only count unique UMIs for each location.\n");
    fprintf(stderr, " -bed      [BED]      Target BED region file. If the strand in column six set, only count reads with the same strand.\n");
    //fprintf(stderr, " -tags     [TAG]      Only count reads with the defined tags.\n");
    fprintf(stderr, " -o        [FILE]     Output depth file. [stdout].\n");
    fprintf(stderr, " -q        [INT]      Minimal map quality to filter. [20]\n");
    fprintf(stderr, " -@        [INT]      Threads to unpack bam. [4]\n");
    fprintf(stderr, "\n\x1b[31m\x1b[1mNotice\x1b[0m :\n");
    fprintf(stderr, " * Requires sorted and indexed BAM as input.\n");
    fprintf(stderr, " * Compares with `samtools depth`, PISA depth considers UMIs and strand of reads.\n");
    return 1;
}
