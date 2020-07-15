
library(VennDiagram)



##### Integration of K4, K27 #####
setwd("data")

# rep1
K4 <- 
  read.csv("K4_rep1_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- 
  read.csv("K27_rep1_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, file="K4K27_rep1_OMO48h_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

# rep2
K4 <- 
  read.csv("K4_rep2_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- 
  read.csv("K27_rep2_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, file="K4K27_rep2_OMO48h_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

# rep3
K4 <- 
  read.csv("K4_rep3_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- 
  read.csv("K27_rep3_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, file="K4K27_rep3_OMO48h_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

# rep4
K4 <- 
  read.csv("K4_rep4_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- 
  read.csv("K27_rep4_OMO48h_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, file="K4K27_rep4_OMO48h_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Preparation for spline #####
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

K4.rep1 <- cbind(K4.K27.rep1[,11:18],K4.K27.rep1[,11:18],K4.K27.rep1[,11:18])
K27.rep1 <- cbind(K4.K27.rep1[,22:29],K4.K27.rep1[,22:29],K4.K27.rep1[,22:29])
K4.rep2 <- cbind(K4.K27.rep2[,11:18],K4.K27.rep2[,11:18],K4.K27.rep2[,11:18])
K27.rep2 <- cbind(K4.K27.rep2[,22:29],K4.K27.rep2[,22:29],K4.K27.rep2[,22:29])
K4.rep3 <- cbind(K4.K27.rep3[,11:18],K4.K27.rep3[,11:18],K4.K27.rep3[,11:18])
K27.rep3 <- cbind(K4.K27.rep3[,22:29],K4.K27.rep3[,22:29],K4.K27.rep3[,22:29])
K4.rep4 <- cbind(K4.K27.rep4[,11:18],K4.K27.rep4[,11:18],K4.K27.rep4[,11:18])
K27.rep4 <- cbind(K4.K27.rep4[,22:29],K4.K27.rep4[,22:29],K4.K27.rep4[,22:29])

K4 <- cbind(K4.rep1,K4.rep2,K4.rep3,K4.rep4)
K27 <- cbind(K27.rep1,K27.rep2,K27.rep3,K27.rep4)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
time <- c(18, 24, 30, 36, 42, 48, 54, 60)
time <- c(time,48+time,96+time)
time <- c(time,time,time,time)
title <- K4.K27.rep1[,1:7]



##### Calculation of max, min, mean, sd, amp #####
maxK4 <- rep(NA,length=nrow(K4))
minK4 <- rep(NA,length=nrow(K4))
meanK4 <- rep(NA,length=nrow(K4))
sdK4 <- rep(NA,length=nrow(K4))
ampK4 <- rep(NA,length=nrow(K4))
maxK27 <- rep(NA,length=nrow(K4))
minK27 <- rep(NA,length=nrow(K4))
meanK27 <- rep(NA,length=nrow(K4))
sdK27 <- rep(NA,length=nrow(K4))
ampK27 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(K4.K27.quadruplicate,maxK4,minK4,meanK4,sdK4,ampK4,
                maxK27,minK27,meanK27,sdK27,ampK27)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(time,K4[i,1:96],spar=0.3)
  spK27 <- smooth.spline(time,K27[i,1:96],spar=0.3)
  length <- 1000
  x <- seq(66,114,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  K4.K27$maxK4[i] <- max(predK4$y)
  K4.K27$minK4[i] <- min(predK4$y)
  K4.K27$meanK4[i] <- mean(predK4$y)
  K4.K27$sdK4[i] <- sd(predK4$y)
  K4.K27$ampK4[i] <- max(predK4$y)-min(predK4$y)
  K4.K27$maxK27[i] <- max(predK27$y)
  K4.K27$minK27[i] <- min(predK27$y)
  K4.K27$meanK27[i] <- mean(predK27$y)
  K4.K27$sdK27[i] <- sd(predK27$y)
  K4.K27$ampK27[i] <- max(predK27$y)-min(predK27$y)
}

write.csv(K4.K27, file="K4K27OMO48h_data_stat.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
K4.K27 <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")

K4.over2 <- subset(K4.K27,maxK4>2)
K27.over2 <- subset(K4.K27,maxK27>2)
K4.K27.over2 <- subset(K4.K27,maxK4>2&maxK27>2)

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.K27.over2, 
          file="K4K27OMO48h_data_stat_maxK4K27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.over2, 
          file="K4K27OMO48h_data_stat_maxK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27OMO48h_data_stat_maxK27over2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27OMO48h_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27OMO48h_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)

K4.K27.over1 <- K4.K27.over1[,1:ncol(K4.over1)]
colnames <- vector(length=ncol(K4.K27.over1))
list <- strsplit(colnames(K4.K27.over1), "\\.x")
for(i in 1:ncol(K4.K27.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.K27.over1) <- colnames

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)

write.csv(K4.over1, 
          file="K4K27OMO48h_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27.over1, 
          file="K4K27OMO48h_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over1, 
          file="K4K27OMO48h_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)



##### Histogram of maximum #####
K4.K27 <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")

maxK4.mod <- rep(NA, length=nrow(K4.K27))
for(i in 1:nrow(K4.K27)){
  if(K4.K27$maxK4[i] <= 8){
    maxK4.mod[i] <- K4.K27$maxK4[i]
  }else{maxK4.mod[i] <- 8.001}
}

maxK27.mod <- rep(NA, length=nrow(K4.K27))
for(i in 1:nrow(K4.K27)){
  if(K4.K27$maxK27[i] <= 8){
    maxK27.mod[i] <- K4.K27$maxK27[i]
  }else{maxK27.mod[i] <- 8.001}
}


pdf("../figs/histogram/histogram_OMO48h_maxK4mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(maxK4.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,8), ylim=c(0,5000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,4000,2000),labels=F)
arrows(-8*0.04,0,-8*0.04,5000,length=0)
mtext(c("0","2","4","6","> 8"), side=1, at=seq(0,8,2),line=-0.3)
mtext(c("0","2,000","4,000"), side=2, at=seq(0,4000,2000),las=1,line=0.15)
mtext(expression(paste(log[2], "(maximum rpkm)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K4me3",side=3,line=0)
text(4.4,3500,"19,571",col="red")
arrows(2,0,2,4000,length=0,col="red")
arrows(2,3500,3,3500,length=0.05,col="red")
dev.off()

pdf("../figs/histogram/histogram_OMO48h_maxK27mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(maxK27.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,8), ylim=c(0,8000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,8000,4000),labels=F)
mtext(c("0","2","4","6","> 8"), side=1, at=seq(0,8,2),line=-0.3)
mtext(c("0","4,000","8,000"), side=2, at=seq(0,8000,4000),las=1,line=0.15)
mtext(expression(paste(log[2], "(maximum rpkm)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K27me3",side=3,line=0)
text(4.4,4000,"9,199",col="red")
arrows(2,0,2,8000,length=0,col="red")
arrows(2,4000,3,4000,length=0.05,col="red")
dev.off()



##### Histogram of amplitude #####
K4 <- read.csv("K4K27OMO48h_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27OMO48h_data_stat_maxK27over2.csv",header=T,sep=",")
nrow(K4)
nrow(K27)

ampK4.mod <- rep(NA, length=nrow(K4))
for(i in 1:nrow(K4)){
  if(K4$ampK4[i] <= 2){
    ampK4.mod[i] <- K4$ampK4[i]
  }else{ampK4.mod[i] <- 2.001}
}

ampK27.mod <- rep(NA, length=nrow(K27))
for(i in 1:nrow(K27)){
  if(K27$ampK27[i] <= 2){
    ampK27.mod[i] <- K27$ampK27[i]
  }else{ampK27.mod[i] <- 2.001}
}

pdf("../figs/histogram/histogram_OMO48h_ampK4mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampK4.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,2), ylim=c(0,6000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,6000,2000),labels=F)
mtext(c("0","1","> 2"), side=1, at=seq(0,2,1),line=-0.3)
mtext(c("0","2,000","4,000","6,000"), side=2, at=seq(0,6000,2000),las=1,line=0.15)
mtext(expression(paste(log[2], "(diel amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K4me3",side=3,line=0)
text(1.05,1400,"1,089",col="red")
arrows(0.5,0,0.5,6000,length=0,col="red")
arrows(0.5,1400,0.8,1400,length=0.05,col="red")
dev.off()

pdf("../figs/histogram/histogram_OMO48h_ampK27mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampK27.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,2), ylim=c(0,3000),breaks=24,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,3000,1000),labels=F)
mtext(c("0","1","> 2"), side=1, at=seq(0,2,1),line=-0.3)
mtext(c("0","1,000","2,000","3,000"), side=2, at=seq(0,3000,1000),las=1,line=0.15)
mtext(expression(paste(log[2], "(diel amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K27me3",side=3,line=0)
text(0.95,700,"44",col="red")
arrows(0.5,0,0.5,3000,length=0,col="red")
arrows(0.5,700,0.8,700,length=0.05,col="red")
dev.off()


# over1
pdf("../figs/histogram/histogram_OMO48h_ampK4mod_over1.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampK4.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,2), ylim=c(0,6000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,6000,2000),labels=F)
mtext(c("0","1","> 2"), side=1, at=seq(0,2,1),line=-0.3)
mtext(c("0","2,000","4,000","6,000"), side=2, at=seq(0,6000,2000),las=1,line=0.15)
mtext(expression(paste(log[2], "(diel amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K4me3",side=3,line=0)
text(1.5,1400,"139",col="red")
arrows(1,0,1,6000,length=0,col="red")
arrows(1,1400,1.2,1400,length=0.05,col="red")
dev.off()

pdf("../figs/histogram/histogram_OMO48h_ampK27mod_over1.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampK27.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,2), ylim=c(0,3000),breaks=24,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,3000,1000),labels=F)
mtext(c("0","1","> 2"), side=1, at=seq(0,2,1),line=-0.3)
mtext(c("0","1,000","2,000","3,000"), side=2, at=seq(0,3000,1000),las=1,line=0.15)
mtext(expression(paste(log[2], "(diel amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K27me3",side=3,line=0)
text(1.5,700,"2",col="red")
arrows(1,0,1,3000,length=0,col="red")
arrows(1,700,1.2,700,length=0.05,col="red")
dev.off()



# Venn diagram
K4 <- read.csv("K4K27OMO48h_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27OMO48h_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:ncol(K4.over1)]

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)


pdf("../figs/edgeR/venn_K4K27_dielampK4K27over1.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4.over1), area2=nrow(K27.over1), cross.area=nrow(K4.K27.over1), 
  cex=1, category=c("ampK4 over1","ampK27 over1"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans", 
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()



##### Ordering of genes according to amplitude of K4 #####
K4 <- read.csv("K4K27OMO48h_data_stat_maxK4over2_ampK4over0.5.csv",header=T,sep=",")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
write.csv(K4.amp.order, 
          file="K4K27OMO48h_data_stat_maxK4over2_ampK4over0.5_ampK4order.csv",
          quote=F, row.names=F)



##### Calculation of mean level #####
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

K4.rep1 <- cbind(K4.K27.rep1[,11:18],K4.K27.rep1[,11:18],K4.K27.rep1[,11:18])
K4.rep2 <- cbind(K4.K27.rep2[,11:18],K4.K27.rep2[,11:18],K4.K27.rep2[,11:18])
K4.rep3 <- cbind(K4.K27.rep3[,11:18],K4.K27.rep3[,11:18],K4.K27.rep3[,11:18])
K4.rep4 <- cbind(K4.K27.rep4[,11:18],K4.K27.rep4[,11:18],K4.K27.rep4[,11:18])

K4 <- cbind(K4.rep1,K4.rep2,K4.rep3,K4.rep4)
K4 <- as.matrix(K4)
time <- c(18, 24, 30, 36, 42, 48, 54, 60)
time <- c(time,48+time,96+time)
time <- c(time,time,time,time)
title <- K4.K27.rep1[,1:7]

K4.mean <- matrix(NA,nrow=nrow(K4),ncol=8)
colnames(K4.mean) <- 
	c("K4_D1T18","K4_D1T24","K4_D2T6","K4_D2T12","K4_D2T18","K4_D2T24","K4_D3T6","K4_D3T12")

K4.max.time.from0.00 <- rep(NA,length=nrow(K4))
K4.min.time.from0.00 <- rep(NA,length=nrow(K4))

K4.res <- cbind(title, K4.mean, 
                K4.max.time.from0.00, K4.min.time.from0.00)

for(i in 1:nrow(K4.res)){
  spK4 <- smooth.spline(time,K4[i,1:96],spar=0.3)
  length <- 1000
  x <- seq(66,114,length=length)
  predK4 <- predict(spK4,x)  
  
  # mean
  K4.res$K4_D1T18[i] <- predK4$y[predK4$x==66]
  K4.res$K4_D1T24[i] <- predK4$y[which.min(abs(predK4$x-72))]
  K4.res$K4_D2T6[i] <- predK4$y[which.min(abs(predK4$x-78))]
  K4.res$K4_D2T12[i] <- predK4$y[which.min(abs(predK4$x-84))]
  K4.res$K4_D2T18[i] <- predK4$y[which.min(abs(predK4$x-90))]
  K4.res$K4_D2T24[i] <- predK4$y[which.min(abs(predK4$x-96))]
  K4.res$K4_D3T6[i] <- predK4$y[which.min(abs(predK4$x-102))]
  K4.res$K4_D3T12[i] <- predK4$y[which.min(abs(predK4$x-108))]


  # Time from 0.00
  if(max(predK4$y)>0){
  	
    if((66+6 <= predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <66+6+24)){
      K4.res$K4.max.time.from0.00[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 66 - 6
    }else if((66+6+24 <= predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <=66+48)){
      K4.res$K4.max.time.from0.00[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 66 - 6 - 24
    }else{
      K4.res$K4.max.time.from0.00[i] <- 
        predK4$x[predK4$y == max(predK4$y)] + 24 - 66 - 6
    }
      
    if((66+6 <= predK4$x[predK4$y == min(predK4$y)]) &&
       (predK4$x[predK4$y == min(predK4$y)] <66+6+24)){
      K4.res$K4.min.time.from0.00[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 66 - 6
    }else if((66+6+24 <= predK4$x[predK4$y == min(predK4$y)]) &&
       (predK4$x[predK4$y == min(predK4$y)] <=66+48)){
      K4.res$K4.min.time.from0.00[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 66 - 6 - 24
    }else{
      K4.res$K4.min.time.from0.00[i] <- 
        predK4$x[predK4$y == min(predK4$y)] + 24 - 66 - 6
    }
  
  }else{}
  
}

write.csv(K4.res, file="K4OMO48h_mean.csv", quote=F, row.names=F)

