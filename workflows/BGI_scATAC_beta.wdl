workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String ref
  String config
  String ?lib
  String tssFile
  String Rscript
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
    ref=ref,
    root=root
  }
  call sortBam {
    input:
    lib=lib,
    bam=fastq2bam.bam,
    root=root,
    outdir=outdir
  }
  call fragIndex {
    input:
    lib=lib,
    root=root,
    fragFile=sortBam.fragFile,
    outdir=outdir
  }
  call summary {
    input:
    lib=lib,
    tssFile=tssFile,
    Rscript=Rscript,
    bam=sortBam.finalBam,    
    root=root,
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
    ${root}/bin/PISA parse -t 10 -q 20 -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -report ${outdir}/report/sequencing_report.csv ${fastq1} ${fastq2} |\
    ${root}/bin/dyncut -1 ${outdir}/temp/reads.fq -l 30 -report ${outdir}/report/trim.csv /dev/stdin
  }
  output {
    String fastq="${outdir}/temp/reads.fq"
    File fastq0="${outdir}/temp/reads.fq"
    File trimReport="{outdir}/report/trim.csv"
  }
}

task fastq2bam {
  String fastq  
  String outdir
  String ref
  String root
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/bin/bwa mem -t 10 -p ${ref} ${fastq} > ${outdir}/temp/aln.sam && \
    ${root}/bin/PISA sam2bam -o ${outdir}/temp/aln.bam -report ${outdir}/report/alignment_report.csv ${outdir}/temp/aln.sam && rm -f ${outdir}/temp/aln.sam
    
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    File bam1="${outdir}/temp/aln.bam"
    File alnReport="{outdir}/report/alignment_report.csv"
  }
}
task sortBam {
  String bam
  String root
  String outdir
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi
    ${root}/bin/sambamba sort -t 10 -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam && \
    ${root}/bin/PISA rmdup -tag CB -t 10 -o ${outdir}/temp/rmdup.bam && \
    ${root}/bin/sambamba index -t 10 ${outdir}/temp/rmdup.bam && \
    ${root}/bin/bap2 bam -i ${outdir}/temp/rmdup.bam -bt CB -r mm10 --mapq 20 -o ${outdir}/outs -c 20 -n sample    
  }
  output {
    String finalBam="${outdir}/outs/sample.bap.bam"
    File finalBam1="${outdir}/outs/sample.bap.bam"
    String fragFile="${outdir}/outs/sample.fragments.gz"
    File fragFile1="${outdir}/outs/sample.fragments.gz"
  }
}

task fragIndex {
  String fragFile
  String root
  String outdir
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/bin/tabix -p bed ${outdir}/outs/sample.fragments.gz
  }
}

task summary {
  String bam
  String root
  String outdir
  String Rscript
  String tssFile
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${root}/bin/bamCoverage -b ${bam} -o ${outdir}/outs/final.bw --numberOfProcessors 10 -of bigwig && \
    ${root}/bin/computeMatrix reference-point -S ${outdir}/out/final.bw -R ${tssFile} -b 3000 -a 3000 --skipZeros -o ${outdir}/temp/final.matrix && \
    ${root}/bin/plotHeatmap -m ${outdir}/temp/final.matrix -out ${outdir}/report/TSS.pdf --colorMap GnBu && rm -f ${outdir}/temp/final.matrix && \
    ${root}/bin/Rscript ${root}/script/scCAT-ATAC-Frasize.R -B ${outdir}/outs/sample.bap.bam -L sample -O ${outdir}/report/Fragsize.pdf
  }
}
