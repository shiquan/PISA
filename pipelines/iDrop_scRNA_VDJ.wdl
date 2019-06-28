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
  String ?runID
  String config
  String cite_config
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
    runID=runID,
    root=root,
  }
  call sortFastq {
    input:
    lib=lib,
    root=root,
    fastq=parseFastq.fastq,
    outdir=outdir,
    config=cite_config,
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
  call report {
    input:
    root = root,
    lib=lib,
    outdir=outdir,
    expectCell = expectCell,
    Rscript=Rscript,
    matrix=countMatrix.matrix,
    ID=ID
  }
}
task sortFastq {
  String root
  String ?lib
  String fastq
  String outdir
  String config
  command {
    source ${lib}
    ${root}/SingleCellTools fsort -tag CB,UY -t 10 -mem -dropN ${fastq} -o ${outdir}/temp/sorted.fq
    ${root}/SingleCellTools parse -config ${config}  ${outdir}/temp/sorted.fq -cbdis /dev/null -repory /dev/null > ${outdir}/temp/select.fq    
  }
}
task report {
  String root
  String Rscript
  String ?lib
  Int ?expectCell
  String ID
  String outdir
  String matrix
  command {
    source ${lib} 
    ${Rscript} -e 'library(rmarkdown);render("${root}/scripts/iDrop_RNAseq.Report.rmd", output_format="html_document", output_file = "${outdir}/outs/${ID}.html", params = list(lib="${ID}",exp="${default=1000 expectCell}"))'
  }
}
task countMatrix {
  String root
  String list
  String outdir
  String anno
  String ?lib
  command {
    source ${lib}

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
    source ${lib}    
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
    source ${lib}
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
    source ${lib}

    ${root}/SingleCellTools parse -t 15 -f -q 20 -dropN -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -run ${default="1" runID} -report ${outdir}/temp/sequencing_report.txt ${fastq1} ${fastq2}  > ${outdir}/temp/reads.fq
    head -n 50000 ${outdir}/temp/barcode_counts_raw.txt |cut -f1 > ${outdir}/temp/barcode_raw_list.txt
  }
  output {
    String count="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq"
    String sequencingReport="${outdir}/temp/sequencing_report.txt"
    String rawlist="${outdir}/temp/barcode_raw_list.txt"
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
    source ${lib}

    ${STAR} --outSAMunmapped Within --outStd SAM --runThreadN 20 --genomeDir ${refdir} --readFilesIn ${fastq} | ${root}/SingleCellTools sam2bam -k -o ${outdir}/temp/aln.bam -report ${outdir}/temp/alignment_report.txt  /dev/stdin
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
    source ${lib}

    ${sambamba} sort -t 20 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam
    ${root}/SingleCellTools anno -gtf ${gtf} -o ${outdir}/temp/anno.bam ${outdir}/temp/sorted.bam
  }
  output {
    String anno="${outdir}/temp/anno.bam"
  }
}
