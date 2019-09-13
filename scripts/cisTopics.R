#!/usr/bin/env Rscript

suppressMessages({
    library(cisTopic)
    library(data.table)
    library(getopt)
})


arg<-matrix(c("input", "i","1","character","input matrix, gzipped",
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
cisObj <- createcisTopicObject(mtx, project.name='mouse liver atac')
rm(mtx)
cisObj <- runModels(cisObj, topic=c(2:30,35,40,45,50), seed=987, nCores=10, burnin = 120, iterations = 150, addModels=FALSE)
save(cisObj, file=opt$output)
