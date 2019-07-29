workflow main {
  String root
  String LFR_fastq1
  String LFR_fastq2
  String fastq1
  String fastq2
  String outdir
  String iDrop_3_config
  String beads_config
  String cell_umi_config
  String refdir
  String hisat
  String sambamba
  String gtf
  call makedir {
    input:
    outdir = outdir
  }
  call parse_fastq_3end {
    input:
    fastq1=fastq1,
    fastq2=fastq2,
    config=iDrop_3_config,
    root = root,
    outdir = makedir.out,
    refdir = refdir,
    hisat = hisat,
    sambamba = sambamba,
    gtf = gtf
  }
  call parse_lfr {
    input:
    fastq1 = LFR_fastq1,
    fastq2 = LFR_fastq2,
    config = beads_config,
    segment_config = cell_umi_config,
    root = root,
    outdir = makedir.out,
    refdir = refdir,
    hisat = hisat,
    sambamba = sambamba,
    gtf = gtf
  }
  call merge {
    input:
    iDrop_bam = parse_fastq_3end.bam,
    LFR_bam = parse_lfr.bam,
    sambamba = sambamba,
    root = root,
    outdir = outdir,    
  }  
}
task makedir {
  String outdir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${outdir}/workflowtime.log
    mkdir -p ${outdir}
    mkdir -p ${outdir}/iDrop
    mkdir -p ${outdir}/scLFR
    mkdir -p ${outdir}/outs
  }
  output {
    String out = "${outdir}"
  }
}
task parse_fastq_3end {
  String fastq1
  String fastq2
  String config
  String root
  String outdir
  String refdir
  String hisat
  String sambamba
  String gtf
  command {
    ${root}/SingleCellTools parse -t 15 -f -q 20 -dropN -config ${config} -cbdis ${outdir}/iDrop/barcode_counts_raw.txt -report ${outdir}/iDrop/sequencing_report.txt ${fastq1} ${fastq2}  > ${outdir}/iDrop/reads.fq 
    ${hisat} -p 20  -x ${refdir} -q ${outdir}/iDrop/reads.fq | ${root}/SingleCellTools sam2bam -o ${outdir}/iDrop/aln.bam /dev/stdin
    ${sambamba} sort -t 20 -o ${outdir}/iDrop/sorted.bam ${outdir}/iDrop/aln.bam
    ${root}/SingleCellTools anno -gtf ${gtf} -o ${outdir}/iDrop/annotated.bam ${outdir}/iDrop/sorted.bam 2>>${outdir}/workflowtime.log
    ${sambamba} index ${outdir}/iDrop/annotated.bam 
    echo "[`date +%F` `date +%T`] iDrop workflow finished" >> ${outdir}/workflowtime.log
  }
  output {
    String bam = "${outdir}/iDrop/annotated.bam"
  }
}
task parse_lfr {
  String fastq1
  String fastq2
  String config
  String segment_config
  String root
  String outdir
  String refdir
  String hisat
  String sambamba
  String gtf
  command {
    ${root}/SingleCellTools parse -t 15 -q 20 -dropN -config ${config} -cbdis ${outdir}/scLFR/beads_barcode_dis.txt -report ${outdir}/scLFR/beads_report.txt ${fastq1} ${fastq2} > ${outdir}/scLFR/read.fq 2>> ${outdir}/workflowtime.log
    ${root}/SingleCellTools fsort -dedup -t 10 -mem -tag LB -p ${outdir}/scLFR/read.fq > ${outdir}/scLFR/sorted.fq 2>> ${outdir}/workflowtime.log
    ${root}/SingleCellTools cleanup -t 10 ${outdir}/scLFR/sorted.fq > ${outdir}/scLFR/cleanup.fq 2>> ${outdir}/workflowtime.log
    ${root}/SingleCellTools unitig -ss GTTCTGCG -t 20 -tag LB ${outdir}/scLFR/cleanup.fq > ${outdir}/scLFR/unitig.fq 2>> ${outdir}/workflowtime.log
    ${root}/SingleCellTools segment2 -config ${segment_config} ${outdir}/scLFR/unitig.fq > ${outdir}/scLFR/parsed_reads.fq 2>> ${outdir}/workflowtime.log
    ${hisat} -p 20  -x ${refdir} -q ${outdir}/scLFR/parsed_reads.fq | ${root}/SingleCellTools sam2bam -o ${outdir}/scLFR/aln.bam /dev/stdin
    ${sambamba} sort -t 20 -o ${outdir}/scLFR/sorted.bam ${outdir}/scLFR/aln.bam
    ${root}/SingleCellTools anno -gtf ${gtf}  -ignore-strand  -splice-consider -t 20 -o ${outdir}/scLFR/annotated.bam ${outdir}/scLFR/sorted.bam 2>>${outdir}/workflowtime.log
    ${sambamba} index ${outdir}/scLFR/annotated.bam 
  }
  output {
    String bam="${outdir}/scLFR/annotated.bam"
  }
}
task merge {
  String iDrop_bam
  String LFR_bam
  String root
  String outdir
  String sambamba
  command {
    ${sambamba} merge -o ${outdir}/outs/merged.bam
    ${root}/SingleCellTools corr -tag UY -tags GN -@ 5 -t 20 -o ${outdir}/outs/UMI_corrected.bam  ${outdir}/outs/merged.bam
    ${root}/SingleCellTools impute -impute CB -tags UY,GN -t 20 -o ${outdir}/outs/imputed_CB.bam ${outdir}/outs/UMI_corrected.bam
    ${root}/SingleCellTools impute -impute CB,UY -tags LB,GN -t 20 -o ${outdir}/outs/final.bam ${outdir}/outs/imputed_CB.bam
  }
}
