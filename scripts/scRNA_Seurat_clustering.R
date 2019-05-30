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
cell <- fread(opt$input, header = T)
cell <- as.data.frame(cell)
names <- as.array(as.matrix(cell[1]))
names = gsub("_", "-", names)
rownames(cell) = names
cell = cell[-1]
names = gsub("-", "_", colnames(cell))
colnames(cell) = names

z <- CreateSeuratObject(counts= cell, min.cells = 3, min.features = 200)
z[["percent.mt"]] <- PercentageFeatureSet(z, pattern = "^mt-")
z <- subset(z, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 3)
z <- NormalizeData(z)
z <-FindVariableFeatures(z, selection.method = "vst", nfeatures = 2000)
})

pdf(paste(opts$output,"/RNA_counts.pdf",sep=""))
VlnPlot(z, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

suppressMessages({
z <- NormalizeData(z)
z <-FindVariableFeatures(z, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(z)
top10 <- head(VariableFeatures(z), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
})
pdf(paste(opts$output,"/top10.pdf",sep=""))
plot2
dev.off()

suppressMessages({
all.genes <- rownames(z)
z <- ScaleData(z, features = all.genes)
z <- RunPCA(z, features = VariableFeatures(object = z))
z <- FindNeighbors(z, dims = 1:10)
z <- FindClusters(z, resolution = 0.5)
z <- RunUMAP(z, dims = 1:10)
})
pdf(paste(opts$output,"/RNA_clustering.pdf",sep=""))
DimPlot(z, reduction = "umap")
dev.off()
