workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String refdir
  String hisat2
  String sambamba
  String gtf
  String config
  String Rscript
  Int ?expectCell
  Int ?forceCell
  call makedir {
    input:
    Dir=outdir
  }
  call parseFastq {
    input:
    lib=lib,
    config=config,
    fastq1=fastq1,
    fastq2=fastq2,
    outdir=makedir.Outdir,
    runID=runID,
    root=root,
  }
  call fastq2bam {
    input:
    lib=lib,
    fastq=parseFastq.fastq,
    outdir=outdir,
    STAR=STAR,
    refdir=refdir,
    root=root
  }
  call sortBam {
    input:
    lib=lib,
    bam=fastq2bam.bam,
    sambamba=sambamba,
    gtf=gtf,
    root=root,
    outdir=outdir
  }
  call cellCount {
    input:
    lib=lib,
    bam=sortBam.anno,
    outdir=outdir,
    root=root,
    rawlist=parseFastq.rawlist
  }    
  call cellCalling {
    input:
    lib=lib,
    root=root,
    count=cellCount.count,
    outdir=outdir,
    Rscript=Rscript,
    expectCell=expectCell,
    forceCell=forceCell,  
  }
  call countMatrix {
    input:
    lib=lib,
    root=root,
    list=cellCalling.list,
    outdir=outdir,
    anno=sortBam.anno
  }

}

task countMatrix {
  String root
  String list
  String outdir
  String anno
  String ?lib
  command {

    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/SingleCellTools count -@ 20 -tag CB -anno_tag GN -umi UY -o ${outdir}/outs/count_mtx.tsv -list ${list} ${anno}
    echo "[`date +%F` `date +%T`] workflow end" >> ${outdir}/workflowtime.log
  }
  output {
    String matrix = "${outdir}/outs/count_mtx.tsv"
  }
}
task cellCalling {
  String count
  String outdir
  String Rscript
  String root
  Int ?expectCell
  Int ?forceCell
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${Rscript} ${root}/scripts/scRNA_cell_calling.R -i ${count} -o ${outdir}/outs -e ${default=1000 expectCell} -f ${default=0 forceCell}
  }
  output {
    String list="${outdir}/outs/cell_barcodes.txt"
  }
}

task cellCount {
  String bam
  String outdir
  String root
  String rawlist
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/SingleCellTools count -tag CB -@ 20 -anno_tag GN -umi UY -o ${outdir}/outs/count_mtx_raw.tsv -count ${outdir}/temp/cell_stat.txt -list ${rawlist} ${bam}
  }
  output {
    String count="${outdir}/temp/cell_stat.txt"
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
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/SingleCellTools parse -t 15 -f -q 20 -dropN -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -report ${outdir}/temp/sequencing_report.txt ${fastq1} ${fastq2}  > ${outdir}/temp/reads.fq

    head -n 50000 ${outdir}/temp/barcode_counts_raw.txt |cut -f1 > ${outdir}/temp/barcode_raw_list.txt

    ${root}/SingleCellTools fsort -tag CB,UY -dedup -list  ${outdir}/temp/barcode_raw_list.txt -dropN -mem -o ${outdir}/temp/sorted.fq ${outdir}/temp/reads.fq

    ${root}/SingleCellTools assem -t 30 -tag CB,UY ${outdir}/temp/sorted.fq > ${outdir}/temp/assem.fq
  }
  output {
    String count="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/assem.fq"
    String sequencingReport="${outdir}/temp/sequencing_report.txt"
    String rawlist="${outdir}/temp/barcode_raw_list.txt"
  }
}

task fastq2bam {
  String fastq  
  String outdir
  String hisat2
  String refdir
  String root
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${STAR}  --outStd SAM --runThreadN 20 --genomeDir ${refdir} --readFilesIn ${fastq} --outFileNamePrefix ${outdir}/temp/ | ${root}/SingleCellTools sam2bam -o ${outdir}/temp/aln.bam -report ${outdir}/temp/alignment_report.txt /dev/stdin
    
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    String alnReport="{outdir}/temp/alignment_report.txt"
  }
}
task sortBam {
  String bam
  String sambamba
  String root
  String outdir
  String gtf
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${sambamba} sort -t 20 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam
    ${root}/SingleCellTools anno -gtf ${gtf} -o ${outdir}/outs/annotated.bam -report ${outdir}/temp/annotated_report.txt ${outdir}/temp/sorted.bam
  }
  output {
    String anno="${outdir}/outs/annotated.bam"
  }
}
