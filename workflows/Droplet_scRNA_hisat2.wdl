workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String refdir
  String hisat2
  String sambamba
  String Rscript
  String ID
  String gtf
  String config
  String ?lib
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
    root=root,
  }
  call fastq2bam {
    input:
    lib=lib,
    fastq=parseFastq.fastq,
    outdir=outdir,
    hisat2=hisat2,
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

    ${root}/PISA count -@ 10 -tag CB -anno_tag GN -umi UB -o ${outdir}/outs/gene_mtx.tsv -list ${list} ${anno}
    gzip -f ${outdir}/outs/gene_mtx.tsv
    echo "[`date +%F` `date +%T`] workflow end" >> ${outdir}/workflowtime.log
  }
  output {
    String matrix = "${outdir}/outs/gene_mtx.tsv"
    File matrix0 = "${outdir}/outs/gene_mtx.tsv"
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

    ${Rscript} ${root}/scripts/scRNA_cell_calling.R -i ${count} -o ${outdir}/outs -e ${default=0 expectCell} -f ${default=0 forceCell}
  }
  output {
    String list="${outdir}/outs/cell_barcodes.txt"
    File list0="${outdir}/outs/cell_barcodes.txt"
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

    ${root}/PISA corr -tag UB -@ 10 -tags-block CB,GN -o ${outdir}/outs/final.bam ${bam} && ${root}/PISA attrcnt -tag CB -tags UB,GN -dedup -list ${rawlist} -@ 10 -o ${outdir}/temp/cell_stat.txt ${outdir}/outs/final.bam
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
    mkdir -p ${Dir}/report
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
  String root
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/PISA parse -t 10 -q 20 -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -report ${outdir}/report/sequencing_report.csv ${fastq1} ${fastq2}  > ${outdir}/temp/reads.fq

    head -n 50000 ${outdir}/temp/barcode_counts_raw.txt |cut -f1 > ${outdir}/temp/barcode_raw_list.txt
  }
  output {
    String count="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq"
    String sequencingReport="${outdir}/report/sequencing_report.csv"
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

    ${hisat2} -x ${refdir} -U ${fastq} -p 10 1> ${outdir}/temp/aln.sam && \
    ${root}/PISA sam2bam -o ${outdir}/temp/aln.bam -report ${outdir}/report/alignment_report.csv ${outdir}/temp/aln.sam && rm -f ${outdir}/temp/aln.sam
    
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    File bam0="${outdir}/temp/aln.bam"
    String alnReport="{outdir}/report/alignment_report.csv"
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

    ${sambamba} sort -t 10 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam
    ${root}/PISA anno -gtf ${gtf} -o ${outdir}/temp/anno.bam -report ${outdir}/report/anno_report.csv ${outdir}/temp/sorted.bam
  }
  output {
    String anno="${outdir}/outs/anno.bam"
  }
}
