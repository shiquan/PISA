#include <stdio.h>


int fastq_parse_usage()
{
    fprintf(stderr, "* Parse cell barcode and UMI string from raw FASTQ.\n");
    fprintf(stderr, "parse [options] lane1_1.fq.gz,lane02_1.fq.gz  lane1_2.fq.gz,lane2_2.fq.gz\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -1       [fastq]   Read 1 output.\n");
    fprintf(stderr, " -2       [fastq]   Read 2 output.\n");
    fprintf(stderr, " -config  [json]    Configure file in JSON format. Required.\n");
    fprintf(stderr, " -run     [string]  Run code, used for different library.\n");
    fprintf(stderr, " -cbdis   [file]    Read count per cell barcode.\n");
    fprintf(stderr, " -p                 Read 1 and read 2 interleaved in the input file.\n");
    // fprintf(stderr, " -t       [INT]     Thread.\n");
    // fprintf(stderr, " -r       [INT]     Records per chunk. [10000]\n");
    // fprintf(stderr, " -dis     [file]    Barcode distribution count.\n");
    fprintf(stderr, " -f                 Filter reads on DNBSEQ standard (2 bases < q10 at first 15 bases).\n");
    fprintf(stderr, " -q       [INT]     Drop reads if average sequencing quality below this value.\n");
    fprintf(stderr, " -dropN             Drop reads if N base in sequence or barcode.\n");
    fprintf(stderr, " -report  [csv]     Summary report.\n");
    fprintf(stderr, "\n");
    return 1;
}
int fsort_usage()
{
    fprintf(stderr, "* Sort reads by tags and deduplicate.\n");
    fprintf(stderr, "fastq-sort [options] in.fq\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, " -tag     [TAGS]     Tags, such as CB,UR. Order of these tags is sensitive.\n");
    fprintf(stderr, " -dedup              Remove dna copies with same tags. Only keep reads have the best quality.\n");
    fprintf(stderr, " -dup-tag [TAG]      Tag name of duplication counts. Use with -dedup only. [DU]\n");
    fprintf(stderr, " -list    [file]     White list for first tag, usually for cell barcodes.\n");
    fprintf(stderr, " -@       [INT]      Threads to compress file.\n");
    fprintf(stderr, " -o       [fq.gz]    bgzipped output fastq file.\n");
    fprintf(stderr, " -m       [mem]      Memory per thread. [1G]\n");
    fprintf(stderr, " -p                  Input fastq is smart pairing.\n");
    fprintf(stderr, " -T       [prefix]   Write temporary files to PREFIX.nnnn.tmp\n");
    fprintf(stderr, " -report  [csv]      Summapry report.\n");
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
    fprintf(stderr, " -dis     [file]     Assembled length distribution.\n");
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
    fprintf(stderr, "* Parse FASTQ+ read name and convert SAM to BAM.\n");
    fprintf(stderr, "sam2bam [options] in.sam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o       [BAM]       Output file [stdout].\n");
    fprintf(stderr, " -mito    [string]    Mitochondria name. Used to stat ratio of mitochondria reads.\n");
    fprintf(stderr, " -maln    [BAM]       Export mitochondria reads into this file instead of standard output file.\n");
    fprintf(stderr, " -@       [INT]       Threads to compress bam file.\n");
    fprintf(stderr, " -report  [csv]       Alignment report.\n");
    // fprintf(stderr, " -skip-sec            Skip secondary alignments in the output.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note :\n");
    fprintf(stderr, "* Reads map to multiple loci usually be marked as low quality and be filter at downstream analysis.\n");
    fprintf(stderr, "  But for RNAseq library, if reads map to a exonic locus but also align to 1 or more non exonic loci,\n");
    fprintf(stderr, "  the exonic locus can be prioritized as primary alignments, and mapping quality adjust to 255. Tag\n");
    fprintf(stderr, "  MM:i:1 will also be added for this record. Following options used to adjust mapping quality.\n");
    fprintf(stderr, "* Input SAM need be sorted by read name, and aligner should output all hits of a read in this SAM.\n");
    fprintf(stderr, " -adjust-mapq         Enable adjust mapping quality score.\n");
    fprintf(stderr, " -gtf     [GTF]       GTF annotation file. This file is required to check the exonic regions.\n");
    fprintf(stderr, " -qual    [255]       Updated quality score.\n");
    fprintf(stderr, "\n");
    return 1;    
}
/*
int rmdup_usage()
{
    fprintf(stderr, "* Deduplicate PCR reads based on read position and barcodes.\n");
    fprintf(stderr, "bam_rmdup [options] in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, "   -tag   [TAGS]       Sample tag, cell barcode tag, and/or UMI tag. RG,CB,UR\n");
    fprintf(stderr, "   -t     [INT]        Threads.\n");
    fprintf(stderr, "   -@     [INT]        Threads to unpack BAM.\n");
    fprintf(stderr, "   -o     [BAM]        Output bam.\n");
    fprintf(stderr, "   -r     [INT]        Records per thread chunk. Default is 10000000.\n");
    fprintf(stderr, "   -k                  Keep duplicates, make flag instead of remove them.\n");
    fprintf(stderr, "\n");
    return 1;
}
*/
int pick_usage()
{
    fprintf(stderr, "* Pick alignment records within barcode list.\n");
    fprintf(stderr, "PickBam [options] in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -list    [file]       Barcode white list.\n");
    fprintf(stderr, " -tag     [TAG]        Barcode tag.\n");
    fprintf(stderr, " -o       [BAM]        Output file.\n");
    fprintf(stderr, " -@       [INT]        Threads to unpack BAM.\n");
    fprintf(stderr, "\n");
    return 1;
}
int anno_usage()
{
    fprintf(stderr, "* Annotate bam records with overlapped function regions. Such as gene, trnascript etc.\n");
    fprintf(stderr, "anno_bam -bed peak.bed -tag PK -o anno.bam in.bam\n");
    fprintf(stderr, "anno_bam -gtf genes.gtf -o anno.bam in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o        [BAM]       Output bam file.\n");
    fprintf(stderr, " -report   [csv]       Summary report.\n");
    fprintf(stderr, " -@        [INT]       Threads to compress bam file.\n");
    fprintf(stderr, " -q        [0]         Map Quality Score cutoff. MapQ smaller and equal to this value will not be annotated.\n");
    fprintf(stderr, "\nOptions for BED file :\n");
    fprintf(stderr, " -bed      [BED]       Function regions. Three or four columns bed file. Col 4 could be empty or names of this region.\n");
    fprintf(stderr, " -tag      [TAG]       Attribute tag name. Set with -bed.\n");
    fprintf(stderr, "\nOptions for mixed samples.\n");
    fprintf(stderr, " -chr-species  [file]  Chromosome name and related species binding list.\n");
    fprintf(stderr, " -btag     [TAG]       Species tag name. Set with -chr-species.\n");
    fprintf(stderr, "\nOptions for GTF file :\n");
    fprintf(stderr, " -gtf      [GTF]       GTF annotation file. gene_id,transcript_id is required for each record.\n");
    fprintf(stderr, " -tags     [TAGS]      Attribute names. Default is TX,AN,GN,GX,RE.\n");
    fprintf(stderr, " -ignore-strand        Ignore strand of transcript in GTF. Reads mapped to antisense transcripts will also be annotated.\n");
    fprintf(stderr, " -splice               Reads covered exon-intron edge will also be annotated.\n");
    fprintf(stderr, " -intron               Reads covered intron regions will be annotated.\n");
    fprintf(stderr, " -t        [INT]       Threads to annotate.\n");
    fprintf(stderr, " -chunk    [INT]       Chunk size per thread.\n");
    fprintf(stderr, "\nNotice :\n");
    fprintf(stderr, " * For GTF mode, this program will set tags in default, you could also reset them by -tags.\n");
    fprintf(stderr, "   TX : Transcript id.\n");
    fprintf(stderr, "   AN : Same with TX but set only if read mapped to antisense strand of transcript.\n");
    fprintf(stderr, "   GN : Gene name.\n");
    fprintf(stderr, "   GX : Gene ID.\n");
    fprintf(stderr, "   RE : Region type, E (Exon), I (Intron), C (Exon and Intron), S (junction reads cover isoforms properly), V (ambiguous reads, junction reads cover isoforms unproperly)\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_corr_usage()
{
    fprintf(stderr, "* Correct error prone barcodes based on frequency.\n");
    fprintf(stderr, "PISA corr [options] in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -o        [BAM]       Output bam.\n");
    fprintf(stderr, " -tag      [TAG]       Tag to correct.\n");
    fprintf(stderr, " -new-tag  [TAG]       Create a new tag for corrected barcodes.\n");
    fprintf(stderr, " -tags-block  [TAGS]   Tags to define read group. For example, if set to GN (gene), reads in the same gene will be grouped together.\n");
    fprintf(stderr, " -cr                   Enable CellRanger like UMI correction method.\n");
    fprintf(stderr, " -@        [INT]       Thread to compress BAM file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Demo : \n");
    fprintf(stderr, " // Two groups of reads have same cell barcode (CB) and gene (GN) but their UMIs (UY) differ by only one base. The UMI of less supported\n");
    fprintf(stderr, " // is corrected to the UMI with higher support. UB save the checked or corrected UMI.\n");
    fprintf(stderr, "* PISA corr -tag UY -new-tag UB -tags-block CB,GN in.bam -o corr.bam \n\n");
    fprintf(stderr, " // Same with above. Besides, if two or more groups of reads have same CB and UB but different GN, the GN with the most supporting reads\n");
    fprintf(stderr, " // is kept for UMI counting, and the other read groups are discarded. In case of a tie for maximal read support, all read groups are\n");
    fprintf(stderr, " // discarded, as the gene cannot be confidently assigned (Cell Ranger method).\n");
    fprintf(stderr, "* PISA corr -cr -tag UY -new-tag UB -tags-block CB,GN in.bam -o corr.bam \n\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_attr_usage()
{
    fprintf(stderr, "* Count the frequency of tag values.\n");
    fprintf(stderr, "PISA count in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -cb       [TAG]      Cell Barcode, or other tag used for each individual.\n");
    fprintf(stderr, " -list     [file]     Cell barcode white list.\n");
    fprintf(stderr, " -tags     [TAGS]     Tags to count.\n");
    fprintf(stderr, " -dedup               Deduplicate the atrributes in each tag.\n");
    fprintf(stderr, " -all-tags            Only records with all tags be count.\n");
    fprintf(stderr, " -group    [TAG]      Group tag, count all tags for each group seperately.\n");
    fprintf(stderr, " -o        [file]     Output count table.\n");
    fprintf(stderr, " -q        [INT]      Map Quality to filter bam.\n");
    fprintf(stderr, " -no-header           Ignore header in the output.\n");
    fprintf(stderr, " -@        [INT]      Thread to unpack bam.\n");
    
    fprintf(stderr, "\n");
    return 1;
}

int bam_extract_usage()
{
    fprintf(stderr, "* Extract tag values from alignments.\n");    
    fprintf(stderr, "bam_extract_tags[options] in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tags     [TAGS]     Tags to be extracted.\n");
    fprintf(stderr, " -o        [file]     Output file. tsv format\n");
    fprintf(stderr, " -n                   Print read name.\n");
    fprintf(stderr, "\n");
    return 1;
}

int bam_count_usage()
{
    fprintf(stderr, "* Count reads or fragments matrix for single-cell datasets.\n");
    fprintf(stderr, "CountMatrix[options] aln.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -tag      [TAG]      Cell barcode tag.\n");
    fprintf(stderr, " -anno-tag [TAG]      Annotation tag, gene or peak.\n");
    fprintf(stderr, " -list     [file]     Barcode white list, used as column names at matrix. If not set, all barcodes will be count.\n");
    //fprintf(stderr, " -o        [file]     Output matrix.\n");
    fprintf(stderr, " -outdir   [DIR]      Output matrix in MEX format into this fold.\n");
    fprintf(stderr, " -umi      [TAG]      UMI tag. Count once if more than one record has same UMI in one gene or peak.\n");
    fprintf(stderr, " -one-hit             Skip if a read hits more than 1 gene or peak.\n");
    fprintf(stderr, " -corr                Enable correct UMIs. Similar UMIs defined as amming distance <= 1.\n");
    fprintf(stderr, " -q        [INT]      Minimal map quality to filter. Default is 20.\n");
    fprintf(stderr, " -@        [INT]      Threads to unpack BAM.\n");
    fprintf(stderr,"\n");
    return 1;
}

int bam2fq_usage()
{
    fprintf(stderr, "* Convert BAM into fastq.\n");
    fprintf(stderr, "bam2fq -tag CB,UY in.bam\n");
    fprintf(stderr, "\nOptions :\n");
    fprintf(stderr, " -filter              Filter this record if not all tags existed.\n");
    fprintf(stderr, " -fa                  Output fasta instead of fastq.\n");
    fprintf(stderr, " -o        [fastq]    Output file.\n");
    fprintf(stderr, " -@        [INT]      Threads to unpack BAM.\n");
    fprintf(stderr, "\n");
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
