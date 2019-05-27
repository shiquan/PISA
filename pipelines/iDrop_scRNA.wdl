workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String refdir
  String STARpath
  String sambambapath
  String ID
  String annoBed
  String ?runID
  String config
  String list
  call makedir {
    input:
    Dir=outdir
  }
  call parseFastq {
    input:
    config=config,
    fastq1=fastq1,
    fastq2=fastq2,
    outdir=makedir.Outdir,
    runID=runID,
    root=root
  }
  call fastq2bam {
    input:
    fastq=parseFastq.fastq,
    outdir=outdir,
    STARpath=STARpath,
    refdir=refdir,
    root=root
  }
  call sortBam {
    input:
    bam=fastq2bam.bam,
    sambambapath=sambambapath,
    annoBed=annoBed,
    root=root,
    list=list,
    outdir=outdir
  }
}

task makedir {
  String Dir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${Dir}/workflowtime.log
    mkdir -p ${Dir}
    mkdir -p ${Dir}/outs
    mkdir -p ${Dir}/temp
  }
  output {
    String Outdir="${Dir}"
  }
}
task parseFastq {
  String config
  String fastq1
  String fastq2
  String outdir
  String ?runID
  String root
  command {
    ${root}/SingleCellTools parse -t 15 -f -q 20 -dropN -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -run ${default="1" runID} -report ${outdir}/temp/sequencing_report.json ${fastq1} ${fastq2} | ${root}/SingleCellTools  trim -report ${outdir}/temp/trim_report.json -mode polyA -p -t 3 > ${outdir}/temp/reads.fq
  }
  output {
    String rawtable="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq"
    String sequencingReport="${outdir}/temp/sequencing_report.json"
  }
}
task fastq2bam {
  String fastq  
  String outdir
  String STARpath
  String refdir
  String root
  command {
    ${STARpath} --outStd SAM --genomeDir ${refdir} --readFilesIn ${fastq} | ${root}/SingleCellTools sam2bam -p MT -o ${outdir}/temp/aln.bam -report ${outdir}/temp/alignment_report.json -maln ${outdir}/temp/mito.bam /dev/stdin
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    String alnReport="{outdir}/temp/alignment_report.json"
  }
}
task sortBam {
  String bam
  String sambambapath
  String root
  String outdir
  String annoBed
  String list
  command {
    ${sambambapath} sort -t 20 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam
    ${root}/SingleCellTools anno -tag GE -bed ${annoBed} -o ${outdir}/temp/anno.bam ${outdir}/temp/sorted.bam
    ${root}/SingleCellTools count -tag CB -anno_tag GE -o ${outdir}/outs/count.mtx -list ${list} ${outdir}/temp/anno.bam
    echo "[`date +%F` `date +%T`] workflow end" >> ${outdir}/workflowtime.log
  }
  output {
    String matrix="${outdir}/outs/count.mtx"
  }
}


