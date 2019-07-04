workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String beads_config
  String runID
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
    counts = parse_fastq.report
  }
  call assem {
    input:
    fastq = sort_fastq.sorted,
    outdir = outdir,
    root = root
  }
}
task assem {
  String fastq
  String outdir
  String root
  command {
    echo "[`date +%F` `date +%T`] assembly start" > ${outdir}/workflowtime.log
    ${root}/SingleCellTools assem -t 20 ${fastq} -tag LB -p > ${outdir}/temp/assem.fq
    echo "[`date +%F` `date +%T`] assembly finished" > ${outdir}/workflowtime.log
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
task parse_fastq {
  String fastq1
  String fastq2
  String outdir
  String root
  String runID
  String beads_config
  command {
    ${root}/SingleCellTools parse -t 15 -q 20 -dropN -config ${beads_config} -cbdis ${outdir}/temp/beads_barcode_dis.txt -run ${default="LFR" runID} -report ${outdir}/temp/beads_report.txt ${fastq1} ${fastq2} > ${outdir}/temp/read.fq
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
  String counts
  command <<<
    awk '{if($2>2) {print $1}}' ${counts}
    ${root}/SingleCellTools fsort -mem -tag LB -o ${outdir}/temp/sorted.fq ${fastq}
  >>>
  output {
    String sorted="${outdir}/temp/sorted.fq"
  }
}
