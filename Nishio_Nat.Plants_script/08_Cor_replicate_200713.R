
library(fields)



##### Correlation between replicates (seasonal ChIP-seq) #####

# Loading data
setwd("data")

K4.K27.rep1 <- read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",header=T,sep=",")
K4.K27.rep2 <- read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",header=T,sep=",")
K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37])

date <- c("6 Nov. 2012","4 Dec. 2012","8 Jan. 2013","5 Feb. 2013","5 Mar. 2013","2 Apr. 2013","30 Apr. 2013","28 May 2013","2 Jul. 2013","30 Jul. 2013","27 Aug. 2013","24 Sep. 2013")


# H3K4me3
pdf("../figs/scatter/scatter_K4_rep1-rep2.pdf",width=7.2,height=2)
set.panel(2,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)

for(i in 1:12){
  par(ps=6)
  plot(K4.rep1[,i], K4.rep2[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K4me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K4me3 rep. 2",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(date[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K4.rep1[,i],K4.rep2[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()


# H3K27me3
pdf("../figs/scatter/scatter_K27_rep1-rep2.pdf",width=7.2,height=2)
set.panel(2,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)

for(i in 1:12){
  par(ps=6)
  plot(K27.rep1[,i], K27.rep2[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9.5),ylim=c(0,9.5), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K27me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K27me3 rep. 2",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(date[i],side=3,line=0)
  arrows(0,0,9.5,9.5,length=0,col="red")
  cor.value <- round(cor(K27.rep1[,i],K27.rep2[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.639,8.972, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()



##### Correlation between replicates (diel ChIP-seq) #####

# Loading data
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_OMO48h_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_OMO48h_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep3 <- 
  read.csv("K4K27_rep3_OMO48h_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep4 <- 
  read.csv("K4K27_rep4_OMO48h_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.quadruplicate <- cbind(K4.K27.rep1,K4.K27.rep2[,8:29],
                              K4.K27.rep3[,8:29],K4.K27.rep4[,8:29])

K4.rep1 <- K4.K27.rep1[,11:18]
K27.rep1 <- K4.K27.rep1[,22:29]
K4.rep2 <- K4.K27.rep2[,11:18]
K27.rep2 <- K4.K27.rep2[,22:29]
K4.rep3 <- K4.K27.rep3[,11:18]
K27.rep3 <- K4.K27.rep3[,22:29]
K4.rep4 <- K4.K27.rep4[,11:18]
K27.rep4 <- K4.K27.rep4[,22:29]

time <- c("18.00 Day1",
          "0.00 Day2","6.00 Day2","12.00 Day2","18.00 Day2",
          "0.00 Day3","6.00 Day3","12.00 Day3")


# H3K4me3
pdf("../figs/scatter/scatter_K4_rep1-rep2_OMO48h.pdf",width=7.2,height=1)
set.panel(1,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)
for(i in 1:8){
  par(ps=6)
  plot(K4.rep1[,i], K4.rep2[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K4me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K4me3 rep. 2",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(time[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K4.rep1[,i],K4.rep2[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()


pdf("../figs/scatter/scatter_K4_rep1-rep3_OMO48h.pdf",width=7.2,height=1)
set.panel(1,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)
for(i in 1:8){
  par(ps=6)
  plot(K4.rep1[,i], K4.rep3[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K4me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K4me3 rep. 3",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(time[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K4.rep1[,i],K4.rep3[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()


pdf("../figs/scatter/scatter_K4_rep1-rep4_OMO48h.pdf",width=7.2,height=1)
set.panel(1,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)
for(i in 1:8){
  par(ps=6)
  plot(K4.rep1[,i], K4.rep4[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K4me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K4me3 rep. 4",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(time[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K4.rep1[,i],K4.rep4[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()


# H3K27me3
pdf("../figs/scatter/scatter_K27_rep1-rep2_OMO48h.pdf",width=7.2,height=1)
set.panel(1,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)
for(i in 1:8){
  par(ps=6)
  plot(K27.rep1[,i], K27.rep2[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K27me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K27me3 rep. 2",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(time[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K27.rep1[,i],K27.rep2[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()


pdf("../figs/scatter/scatter_K27_rep1-rep3_OMO48h.pdf",width=7.2,height=1)
set.panel(1,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)
for(i in 1:8){
  par(ps=6)
  plot(K27.rep1[,i], K27.rep3[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K27me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K27me3 rep. 3",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(time[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K27.rep1[,i],K27.rep3[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()


pdf("../figs/scatter/scatter_K27_rep1-rep4_OMO48h.pdf",width=7.2,height=1)
set.panel(1,8,relax=T)
par(mgp=c(0,0.3,0))
par(mar=c(1.2,1.1,0.7,0.2))
par(xpd=T)
par(cex=1)
for(i in 1:8){
  par(ps=6)
  plot(K27.rep1[,i], K27.rep4[,i], pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,9),ylim=c(0,9), main="",axes=F)
  axis(1, tcl=-0.1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=1,line=-0.35)
  mtext("H3K27me3 rep. 1",side=1,line=0.1)
  #mtext("amplitude",side=1,line=1.1)
  axis(2, tcl=-0.1, las=1, labels=F)
  mtext(c("0","2","4","6","8"),at=seq(0,8,2),side=2,line=0.15,las=1)
  mtext("H3K27me3 rep. 4",side=2,line=0.4)
  #mtext("amplitude",side=2,line=0.8)
  mtext(time[i],side=3,line=0)
  arrows(0,0,9,9,length=0,col="red")
  cor.value <- round(cor(K27.rep1[,i],K27.rep4[,i], method="spearman"), digits=2)
  par(ps=5)
  text(2.5,8.5, paste("r = ",cor.value, sep=""))
  box()
}
set.panel()
dev.off()
