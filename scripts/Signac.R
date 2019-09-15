#!/usr/bin/env Rscript

suppressMessages({
    library(Signac)
    library(Seurat)
    library(data.table)
    library(getopt)
})


arg<-matrix(c("input", "i","1","character","input matrix, gzipped",
              "frag", "f", "1", "character", "input fragment file. Tabix indexed.",
              "output","o","2","character","output R object",
              "help"  ,"h","2","logical",  "This information"),
            byrow=T,ncol=5
            )

opt = getopt(arg)

if( !is.null(opt$help) || is.null(opt$input)){
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

mtx <- as.data.frame(fread(paste("gzip -dc ",opt$input, sep=""),header=T))
rownames(mtx) = mtx$ID
mtx = mtx[,-1]
obj <- CreateSeuratObject(counts = mtx, assay = 'peaks', project = 'ATAC', min.cells = 10)
rm(mtx)


obj <- SetFragments(
  object = obj,
  file = opt$frag
)

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(
  object = obj,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 1:30)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 1:30)
obj <- FindClusters(object = obj, verbose = FALSE)


save(obj, file=opt$output)
