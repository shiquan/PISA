#!/usr/bin/env Rscropt

suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
    library(cowplot)
#    library(Cairo)
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

if (is.null(opt$output)) {
    opt$output<-getwd()
}

bc <- fread(opt$input,header=TRUE)
bc <- as.data.frame(bc)
bc <- subset(bc, bc$nUMI>0)
bc <- bc[order(bc$nUMI, decreasing=T),]
len <- nrow(bc)
                                        #sor = sort(bc$nUMI, decreasing=T)
sor = bc$nUMI
a = log10(1:len)
b = log10(sor)
expect <- 0
cutoff <- 0
m <- 0
if (!is.null(opt$expect)) {
    expect=as.numeric(opt$expect)
    lo <- loess(b~a,span = 0.004,degree = 2)
    c = log10(expect*1.7)
    if(10^c>len){c=log10(expect)}    
    xl <- seq(2,c, (c - 2)/10)
    out = predict(lo,xl)
    infl <- c(FALSE,abs(diff(out)/((c - 2)/10) - -1) == min(abs(diff(out)/((c - 2)/10)- -1)))
    m = 10 ^ out[infl] + 0.5
    m = round(m , digits =0 )
    cutoff<-length(which(sor>=m))
}

if (!is.null(opt$force)) {
    force = as.numeric(opt$force)
    if (force > 0) {
        expect = force
        cutoff = expect
        m = sor[cutoff]
    }
}

tmp<-data.frame(x=1:len,y=sor,cell=c(rep("true",cutoff),rep("noise",len-cutoff)))

cc <- bc[1:cutoff,]
called = sum(cc$nUMI)
in_cell <- round(called/sum(bc$nUMI),digits=4)

small<-bc[which(bc$nUMI>=m),]
CellNum = nrow(small)
UMI_mean = mean(small$nUMI)
UMI_median = median(small$nUMI)
Gene_mean = mean(small$nGene)
Gene_median = median(small$nGene)
UMI_Cell = sum(small$nUMI)

write.table(cc$CELL_BARCODE, file=paste(opt$output,"/cell_barcodes.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

                                        #png(file=paste(opt$output,"/cell_count_summary.png",sep=""), width=700,height=300,res=100,pointsize=10)
png(file=paste(opt$output,"/cell_count_summary.png",sep=""), width=1200,height=400,res=80)
p1 = ggplot(tmp,aes(x=x,y=y)) + xlim(0,10) + ylim(0,9)
p1 = p1 +annotate("text",x=0.2,y=9,label="Estimated Number of Cell:",size=10,hjust=0)
p1 = p1 +annotate("text",x=3,y=6.5,label=CellNum,size=20,hjust=0, colour="red")
p1 = p1 +annotate("text",x=0.2,y=3.5,label="Reads in cell:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=3.5,label=round(in_cell,3),size=6,hjust=1)
p1 = p1 +annotate("segment",x = 0, xend = 10, y = 3, yend = 3,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=2.5,label="Mean UMI counts per cell:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=2.5,label=round(UMI_mean),size=6,hjust=1)
p1 = p1 +annotate("segment", x = 0, xend = 10, y = 2, yend = 2,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=1.5,label="Mean Genes per cell:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=1.5,label=round(Gene_mean),size=6,hjust=1)
p1 = p1 +annotate("segment",x  = 0, xend  = 10, y  = 1, yend = 1,colour = "blue")
p1 = p1 +theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major=element_blank(),panel.grid.minor = element_blank())

p = ggplot(tmp,aes(x=x,y=y))
p = p + geom_line(aes(color=cell),size=2) +scale_color_manual(values=c("#999999","blue"))
p = p + scale_x_log10(name="Barcodes",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + scale_y_log10(name="nUMI",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + theme_bw() + geom_vline(xintercept =cutoff)
#p = p + geom_text(aes(x=10,y=10,label = paste("cell=",cutoff)), color = 'blue',size=6)
#p = p + geom_text(aes(x=10,y=4,label = paste("nUMI=",m)), color = 'blue',size=6)
#p = p + geom_text(aes(x=10,y=1,label = paste("nUIM%=",in_cell)), color = 'blue',size=6)
p = p + theme(legend.position = "none")

                                        #p1 <- ggplot(bc) + geom_boxplot(aes(x=5,y=nUMI), outlier.shape = 8, width=10) + theme_classic(
#p3 <- p3 + geom_jitter(aes(x=sample(1:10,nrow(bc),replace = T),y=nUMI),alpha=0.2,color="blue") #+ scale_y_log10()
p3 <- ggplot(cc) + geom_violin(aes(x=5,y=nUMI),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=nUMI),alpha=0.2,color="blue")
p3 <- p3 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p2 <- ggplot(bc) + geom_boxplot(aes(x=5,y=nGene), outlier.shape = 8, width=10) + theme_classic()
                                        #p2 <- p2 + geom_jitter(aes(x=sample(1:10,nrow(bc),replace = T),y=nGene),alpha=0.2,color="blue") #+ scale_y_log10()
p2 <- ggplot(cc) + geom_violin(aes(x=5,y=nGene),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=nGene),alpha=0.2,color="blue")
p2 <- p2 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

plot_grid(p1, p, p3,p2, rel_widths = c(4,3,1.5,1.5),ncol=4)

dev.off()
