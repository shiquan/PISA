workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String beads_config
  String seg_config
  String runID
  String hisat2
  String refdir
  String sambamba
  String bwa
  call makedir {
    input:
    outdir = outdir
  }
  call parse_fastq {
    input:
    fastq1 = fastq1,
    fastq2 = fastq2,
    outdir = makedir.dir,
    beads_config = beads_config,
    root = root,
    runID = runID
  }
  call sort_fastq {
    input:
    fastq = parse_fastq.fastq,
    outdir = outdir,
    root = root,
    bwa = bwa,
  }
  call assem {
    input:
    fastq = sort_fastq.sorted,
    outdir = outdir,
    root = root
  }
  call align1 {
    input:
    outdir = outdir,
    fastq = assem.assem,
    root = root,
    hisat2 = hisat2,
    refdir = refdir,
    sambamba = sambamba,
    seg_config = seg_config,
  }
  call align2 {
    input:
    outdir = outdir,
    fastq = sort_fastq.sorted,
    root = root,
    hisat2 = hisat2,
    refdir = refdir,
    sambamba = sambamba,
    seg_config = seg_config,
  }
}
task align1 {
  String fastq
  String root
  String hisat2
  String refdir
  String sambamba
  String seg_config
  String outdir
  command {
    ${root}/SingleCellTools segment2 -t 5 -config ${seg_config} ${fastq} > ${outdir}/temp/assem_format.fa
    ${hisat2} -x ${refdir} -f ${outdir}/temp/assem_format.fa | ${root}/SingleCellTools sam2bam -o ${outdir}/temp/assem_aln.bam
    ${sambamba} sort -o ${outdir}/outs/assem.bam
  }
}
task align2 {
  String fastq
  String root
  String hisat2
  String refdir
  String sambamba
  String seg_config
  String outdir
  command {
    ${root}/SingleCellTools segment2 -t 5 -config ${seg_config} ${fastq} > ${outdir}/temp/pemerge_format.fa
    ${hisat2} -x ${refdir} -f ${outdir}/temp/pemerge_format.fa | ${root}/SingleCellTools sam2bam -o ${outdir}/temp/pemerge_aln.bam
    ${sambamba} sort -o ${outdir}/outs/pemerge.bam
  }
}

task assem {
  String fastq
  String outdir
  String root
  command {
    echo "[`date +%F` `date +%T`] assembly start" >> ${outdir}/workflowtime.log
    ${root}/SingleCellTools assem -t 20 ${fastq} -tag LB -p > ${outdir}/temp/assem.fq
    echo "[`date +%F` `date +%T`] assembly finished" >> ${outdir}/workflowtime.log
  }
  output {
    String assem = "${outdir}/temp/assem.fq"
    }
}
task makedir {
  String outdir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${outdir}/workflowtime.log
    mkdir -p ${outdir}
    mkdir -p ${outdir}/outs
    mkdir -p ${outdir}/temp
  }
  output {
    String dir="${outdir}"
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
    ${root}/SingleCellTools parse -t 10 -q 20 -dropN -config ${beads_config} -cbdis ${outdir}/temp/beads_barcode_dis.txt -run ${default="LFR" runID} -report ${outdir}/temp/beads_report.txt ${fastq1} ${fastq2} > ${outdir}/temp/read.fq
  }
  output {
    String count="${outdir}/temp/beads_barcode_dis.txt"
    String fastq="${outdir}/temp/read.fq"
    String report="${outdir}/temp/beads_report.txt"
  }
}
task sort_fastq {
  String fastq
  String outdir
  String root
  String bwa
  command <<<
    ${root}/SingleCellTools fsort -dedup -t 10 -mem -tag LB -p ${fastq} > ${outdir}/temp/sorted.fq
    ${root}/SingleCellTools trim -t 10 -mode Tn5 -d -tail 10 ${outdir}/temp/sorted.fq > ${outdir}/temp/trimmed.fq
  >>>
  output {
    String sorted="${outdir}/temp/sorted.fq"
  }
}
