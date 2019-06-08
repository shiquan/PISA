workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String refdir
  String bwapath
  String sambambapath
  String ID
  String ?runID
  String ?gsize
  String config
  String macspath
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
    bwapath=bwapath,
    refdir=refdir,
    root=root
  }
  call sortBam {
    input:
    bam=fastq2bam.bam,
    sambambapath=sambambapath,
    root=root,
    outdir=outdir
  }
  call callPeak {
    input:
    macspath=macspath,
    bam=sortBam.sorted,
    outdir=outdir,
    root=root,
    ID=ID,
    gsize=gsize
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
    ${root}/SingleCellTools parse -t 15 -f -q 20 -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -run ${default="1" runID} -report ${outdir}/temp/sequencing_report.json ${fastq1} ${fastq2} | ${root}/SingleCellTools  trim -report ${outdir}/temp/trim_report.json -mode Tn5 -p -t 3 > ${outdir}/temp/reads.fq
  }
  output {
    String rawtable="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq"
    String sequencingReport="${outdir}/temp/sequencing_report.txt"
  }
}
task fastq2bam {
  String fastq  
  String outdir
  String bwapath
  String refdir
  String root
  command {
    ${bwapath} mem -t 20 -p ${refdir} ${fastq} | ${root}/SingleCellTools sam2bam -p -o ${outdir}/temp/aln.bam -report ${outdir}/temp/alignment_report.txt -maln ${outdir}/temp/mito.bam /dev/stdin
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    String alnReport="{outdir}/temp/alignment_report.txt"
  }
}
task sortBam {
  String bam
  String sambambapath
  String root
  String outdir
  command {
    ${sambambapath} sort -t 20 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam 
    ${root}/SingleCellTools rmdup -tag CB -t 5 -o ${outdir}/temp/rmdup.bam ${outdir}/temp/sorted.bam
  }
  output {
    String sorted="${outdir}/temp/rmdup.bam"
  }
}
task callPeak {
  String bam
  String root
  String outdir
  String macspath
  String ID
  String ?gsize
  command <<<
    ${macspath} callpeak -t ${bam} -f BAM --keep-dup all --nomodel --shift -100 --extsize 200 -g ${default="mm" gsize} -n ${ID} --outdir ${outdir}/outs
    cut -f1,2,3 ${outdir}/outs/${ID}_peaks.narrowPeak > ${outdir}/temp/peak.bed
    ${root}/SingleCellTools anno -bed ${outdir}/temp/peak.bed -tag PK -o ${outdir}/outs/processed.bam ${bam}
    ${root}/SingleCellTools attrcnt -cb CB -tag PK -o ${outdir}/temp/readcount.report.txt ${outdir}/outs/processed.bam
    awk '{if($1 !~ /CELL_BARCODE/ && $2>1000 && $3/$2>0.1){print $1;}}' ${outdir}/temp/readcount.report.txt > ${outdir}/temp/barcodes_called.txt
    awk '{if($1 !~ /CELL_BARCODE/ && $2>1000 && $3/$2>0.1){print $0;}}' ${outdir}/temp/readcount.report.txt |awk 'BEGIN{rw=0;tar=0;i=0}{rw+=$2; tar+=$3; i++} END{print rw/i,tar/i}' >> ${outdir}/workflowtime.log
    ${root}/SingleCellTools count -tag CB -anno_tag PK -list ${outdir}/temp/barcodes_called.txt -o ${outdir}/outs/count_matrix.txt ${outdir}/outs/processed.bam
    echo "[`date +%F` `date +%T`] workflow end" >> ${outdir}/workflowtime.log
  >>>
}
