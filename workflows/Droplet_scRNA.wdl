workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String refdir
  String STAR
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
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/PISA corr -tag UB -@ 10 -tags-block CB,GN -o ${outdir}/outs/final.bam ${bam} && \
    ${root}/PISA attrcnt -tag CB -tags UB,GN -dedup -@ 10 -o ${outdir}/temp/cell_stat.txt -all-tags ${outdir}/outs/final.bam
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

    # disable multi-threads support
    ${root}/PISA parse -q 20 -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -report ${outdir}/report/sequencing_report.csv ${fastq1} ${fastq2} -1 ${outdir}/temp/reads.fq

  }
  output {
    String count="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq"
    String sequencingReport="${outdir}/report/sequencing_report.csv"
      }
}

task fastq2bam {
  String fastq  
  String outdir
  String STAR
  String refdir
  String root
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${STAR}  --outStd SAM --runThreadN 10 --genomeDir ${refdir} --readFilesIn ${fastq} --outFileNamePrefix ${outdir}/temp/ 1> ${outdir}/temp/aln.sam && \
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
    String anno="${outdir}/temp/anno.bam"
  }
}
