workflow main {
  String root
  String fastq1
  String fastq2
  String refdir
  String STAR
  String sambamba
  String gtf
  String outdir
  String beads_config
  String cell_umi_config
  String umi_config
  String runID
  call makedir {
    input:
    outdir = outdir
  }
  call parse_fastq {
    input:
    fastq1=fastq1,
    fastq2=fastq2,
    outdir=makedir.dir,
    beads_config = beads_config,
    root=root,
    runID=runID
  }
  call parse_cell_umi {
    input:
    root = root,
    read1 = parse_fastq.read1,
    read2 = parse_fastq.read2,
    config= cell_umi_config,
    outdir = outdir,
  }
  call parse_umi {
    input:
    root = root,
    read1 = parse_fastq.read1,
    read2 = parse_fastq.read2,
    config = umi_config,
    outdir = outdir,
  }
  call aln_all_reads {
    input:
    read1 = parse_fastq.read1,
    read2 = parse_fastq.read2,
    gtf = gtf,
    sambamba = sambamba,
    STAR = STAR,
    refdir = refdir,
    outdir = outdir,
    root = root,
  }

  call report {
    input:
    parse_report = parse_fastq.report,
    cell_umi = parse_cell_umi.report,
    umi = parse_umi.report,
    outdir = outdir
  }
}
task report {
  String parse_report
  String cell_umi
  String umi
  String outdir
  command <<<
    cat ${parse_report} > ${outdir}/outs/summary_report.txt
    cat ${cell_umi} >> ${outdir}/outs/summary_report.txt
    cat ${umi}  >> ${outdir}/outs/summary_report.txt
    echo "[`date +%F` `date +%T`] workflow end" > ${outdir}/workflowtime.log
  >>>
}

task aln_all_reads {
  String read1
  String read2
  String STAR
  String gtf
  String sambamba
  String root
  String outdir
  String refdir
  command {
    ${STAR} --outStd SAM --runThreadN 20 --genomeDir ${refdir} --readFilesIn ${read1} ${read2} | ${root}/SingleCellTools sam2bam -o ${outdir}/temp/aln.bam -report ${outdir}/temp/alignment_report.txt -maln ${outdir}/temp/mito.bam /dev/stdin
    ${sambamba} sort -t 20 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam
    ${root}/SingleCellTools anno -gtf ${gtf} -ignore-strand -o ${outdir}/temp/annotated.bam ${outdir}/temp/sorted.bam    
  }
  output {
    String report="${outdir}/temp/alignment_report.txt"
  }
}
task parse_umi {
  String read1
  String read2
  String outdir
  String root
  String config
  command <<<
    ${root}/SingleCellTools segment -config ${config} -1 ${read1} -2 ${read2} -cb LB -pair ${outdir}/temp/umi/UMI_pair.fq -o ${outdir}/temp/umi/segments_UMI.txt
    awk '{if($1 !~ /RNAME/){printf("%s\t%s\n",$1,$3)}}' ${outdir}/temp/umi/segments_UMI.txt | sed -r "s/.*LB\|\|\|(.*)/\1/" | sort -k1,1 -k2,2  |uniq > ${outdir}/temp/umi/beads_umi_uniq.txt
    n=`wc -l ${outdir}/temp/umi/beads_umi_uniq.txt`
    umi_per_beads=`cut -f1 ${outdir}/temp/umi/beads_umi_uniq.txt|uniq -c|awk 'BEGIN{i=0;j=0}{i++; j+=$1}END{print i/j}`
    echo "UMI per Beads : $umi_per_beads" > ${outdir}/temp/cell_barcode_umi/report.txt
    echo "Beads with partly Cell Barcode + UMI : $n" >> ${outdir}/temp/cell_barcode_umi/report.txt
  >>>
  output {
    String report="${outdir}/temp/umi/report.txt"
    String fastq="${outdir}/temp/umi/UMI_pair.fq"
  }
}

task parse_cell_umi {
  String read1
  String read2
  String outdir
  String root
  String config
  command <<<
    ${root}/SingleCellTools segment -config ${config} -1 ${read1} -2 ${read2} -cb LB -pair ${outdir}/temp/cell_barcode_umi/Cell_Barcode_UMI_pair.fq -o ${outdir}/temp/cell_barcode_umi/segments_Cell_Barcode_UMI.txt
    awk '{if($1 !~ /RNAME/){printf("%s\t%s%s\t%s\n",$1,$2,$3,$4)}}' ${outdir}/temp/cell_barcode_umi/segments_Cell_Barcode_UMI.txt | sed -r "s/.*LB\|\|\|(.*)/\1/" | sort -k1,1 -k2,2 -k3,3 |uniq > ${outdir}/temp/cell_barcode_umi/beads_cell_umi_uniq.txt
    cut -f2 ${outdir}/temp/cell_barcode_umi/beads_cell_umi_uniq.txt |sort |uniq -c | awk '{printf("%s\t%d\n",$2,$1)}' > ${outdir}/temp/cell_barcode_umi/cell_barcodes_count.txt
    bc_per_beads=`cut -f1,2 ${outdir}/temp/cell_barcode_umi/beads_cell_umi_uniq.txt |uniq| cut -f1|uniq -c | awk 'BEGIN{i=0;j=0}{i++; j+=$1}END{print i/j}`
    n=`wc -l ${outdir}/temp/cell_barcode_umi/beads_cell_umi_uniq.txt`
    echo "Beads with Cell Barcode + UMI : $n" >> ${outdir}/temp/cell_barcode_umi/report.txt
    echo "Cell Barcode per Beads : $bc_per_beads" > ${outdir}/temp/cell_barcode_umi/report.txt        
  >>>
  output {
    String report="${outdir}/temp/cell_barcode_umi/report.txt"
    String fastq="${outdir}/temp/cell_barcode_umi/Cell_Barcode_UMI_pair.fq"
  }
}
task parse_fastq {
  String fastq1
  String fastq2
  String outdir
  String root
  String runID
  String beads_config
  command {
    ${root}/SingleCellTools parse -t 15 -q 20 -dropN -config ${beads_config} -cbdis ${outdir}/temp/beads_barcode_dis.txt -run ${default="LFR" runID} -report ${outdir}/temp/beads_report.txt ${fastq1} ${fastq2} -1 ${outdir}/temp/read_1.fq -2 ${outdir}/temp/read_2.fq    
  }
  output {
    String count="${outdir}/temp/beads_barcode_dis.txt"
    String read1="${outdir}/temp/read_1.fq"
    String read2="${outdir}/temp/read_2.fq"
    String report="${outdir}/temp/beads_report.txt"
  }
}
task makedir {
  String outdir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${outdir}/workflowtime.log
    mkdir -p ${outdir}
    mkdir -p ${outdir}/outs
    mkdir -p ${outdir}/temp
    mkdir -p ${outdir}/temp/cell_barcode_umi
    mkdir -p ${outdir}/temp/umi
  }
  output {
    String dir="${outdir}"
  }
}
