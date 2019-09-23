#!/usr/bin/env Rscropt

suppressMessages({
    library(ggplot2)
    library(data.table)
    library(Seurat)
    library(getopt)
})

arg<-matrix(c("input", "i","1","character","Path of input directory",
              "output","o","2","character","Path of output",
              "help"  ,"h","2","logical",  "This information"),
            byrow=T,ncol=5
            )

opt = getopt(arg)

if( !is.null(opt$help) || is.null(opt$input) || is.null(opt$expect)  ){
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

ifelse(is.null(opt$output),opt$output<-getwd(),"")

suppressMessages({
mtx <- as.data.frame(fread(paste("gzip -dc ",opt$input, sep=""),header=T))
rownames(mtx) = mtx$ID
mtx = mtx[,-1]

z <- CreateSeuratObject(counts= mtx, min.cells = 3, min.features = 200)
rm(mtx)
z[["percent.mt"]] <- PercentageFeatureSet(z, pattern = "^mt-")
z <- subset(z, subset = nFeature_RNA > 500 & percent.mt < 10)
z <- NormalizeData(z)
z <-FindVariableFeatures(z, selection.method = "vst", nfeatures = 2000)
})

pdf(paste(opts$output,"/RNA_counts.pdf",sep=""))
VlnPlot(z, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

suppressMessages({
z <- NormalizeData(z)
z <-FindVariableFeatures(z, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(z)
z <- ScaleData(z, features = all.genes)
z <- RunPCA(z, features = VariableFeatures(object = z))
z <- FindNeighbors(z, dims = 1:10)
z <- FindClusters(z, resolution = 0.5)
z <- RunUMAP(z, dims = 1:10)
})

cluster_ID=as.data.frame(Idents(object = z))
cluster_cor= as.data.frame(Embeddings(object = z,reduction = "umap"))

coor=cbind(cluster_ID,cluster_cor,z[['percent.mt']],z[['nCount_RNA']],z[['nFeature_RNA']])
colnames(coor) = c("Cluster","UMAP_1","UMAP_2","percent.mt","nUMI","nGene")

write.table(coor, file=paste(opt$output,"/cluster.csv",sep=""), sep=",",col.names=T)
