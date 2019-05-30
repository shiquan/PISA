#!/usr/bin/env Rscropt

suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
})

arg<-matrix(c("input", "i","1","character","Path of input directory",
              "output","o","2","character","Path of output",
              "force","f","1","numeric","Force cells number to analysis. You must set -f or -e to process barcode list.",
              "expect","e","1","numeric","Number of expect cells",
              "help"  ,"h","2","logical",  "This information"),
            byrow=T,ncol=5
            )

opt = getopt(arg)

if( !is.null(opt$help) || is.null(opt$input)){
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

if (is.null(opt$expect) && is.null(opt$force)) {
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

if(is.null(opt$output)) {
    opt$output<-getwd()
}
expect <- 0
if (!is.null(opt$expect)) {
    expect=as.numeric(opt$expect)
}
# force = expect
if (!is.null(opt$force)) {
    expect = as.numeric(opt$force)
}

bc <- fread(opt$input,header=FALSE)
bc <- as.data.frame(bc)
colnames(bc) <- c("bc","counts")
len <- nrow(bc)
a = log10(1:len)
b =log10(as.numeric(bc$counts))

lo <- loess(b~a,span = 0.004,degree = 2)
#expect <- 1000
c = log10(expect*1.7)
if(10^c>len){c=log10(expect)}
# print(c)
xl <- seq(2,c, (c - 2)/10)
out = predict(lo,xl)
infl <- c(FALSE,abs(diff(out)/((c - 2)/10) - -1) == min(abs(diff(out)/((c - 2)/10)- -1)))
m = 10 ^ out[infl] + 0.5
m = round(m , digits =0 )

cutoff<-length(which(bc$counts>=m))

tmp<-data.frame(x=1:len,y=bc$counts,g=c(rep("true",cutoff),rep("noise",len-cutoff)))
cbs <- bc$bc[1:len]
write.table(cbs, file=paste(opt$output,"/cell_barcodes.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
pdf(paste(opt$output,"/cell_calling.pdf",sep=""))
p = ggplot(tmp,aes(x=x,y=y))
p = p +geom_line(aes(color=g),size=2) +scale_color_manual(values=c("#999999","blue"))
p = p +scale_x_log10(name="Barcodes")
p = p +scale_y_log10(name="Raw counts",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + theme_bw() + geom_vline(xintercept =cutoff)
p = p +  geom_text(aes(x=10,y=1,label = paste("cell=",cutoff)), color = 'blue',size=4)
p = p +  geom_text(aes(x=10,y=2,label = paste("counts=",m)), color = 'blue',size=4)
p
dev.off()
