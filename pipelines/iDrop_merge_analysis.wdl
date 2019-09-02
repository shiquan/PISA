workflow main {
  String root
  String bams
  String outdir
  String sambamba
  call makedir {
    input:
    dir=outdir
  }
  call mergebams {
    input:
    bams = bams,    
    outdir = makedir.outdir,
    sambamba = sambamba
  }
  call barcodes {
    input:
    outdir = outdir,
    root = root,
    bam = mergebams.merged
  }
  call count_raw {
    input:
    outdir = outdir,
    root = root,
    bam = mergebams.merged,
    list = barcodes.list
  }
  call cell_calling {
    input:
    outdir = outdir,
    root = root,
    count = count_raw.count,
  }
  call filter_mtx {
    outdir = outdir,
    root = root,
    list = cell_calling.list
  }
}
task mergebams {
  String bams

  
}
task makedir {
  String dir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${Dir}/workflowtime.log
    mkdir -p ${dir}
    mkdir -p ${dir}/outs
    mkdir -p ${dir}/temp
  }
  output {
    String outdir="${dir}"
  }
}
