
library(fields)
library(LSD)
library(calibrate)
library(cocor)

source("functions/ng.Colors.R")



###### Integration of H3K4me3 and H3K27me3 data #####
setwd("data")

# rep1
K4 <- read.csv("K4_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_rep1_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_rep1_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

#rep2
K4 <- read.csv("K4_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_rep2_rpkm+1_log2_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_rep2_rpkm+1_log2_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22],K4.K27.rep1[,11:22],K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37],K4.K27.rep1[,26:37],K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22],K4.K27.rep2[,11:22],K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37],K4.K27.rep2[,26:37],K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27.duplicate[,1:7]



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

K4.K27 <- cbind(K4.K27.duplicate,maxK4,minK4,meanK4,sdK4,ampK4,
                maxK27,minK27,meanK27,sdK27,ampK27)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
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
write.csv(K4.K27, file="K4K27_data_stat.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
K4.K27 <- read.csv("K4K27_data_stat.csv",header=T,sep=",")

K4.over2 <- subset(K4.K27,maxK4>2)
K27.over2 <- subset(K4.K27,maxK27>2)
K4.K27.over2 <- subset(K4.K27,maxK4>2&maxK27>2)

nrow(K4.over2)
nrow(K27.over2)
nrow(K4.K27.over2)

write.csv(K4.K27.over2, 
          file="K4K27_data_stat_maxK4K27over2.csv", 
          quote=F, row.names=F)
write.csv(K4.over2, 
          file="K4K27_data_stat_maxK4over2.csv", 
          quote=F, row.names=F)
write.csv(K27.over2, 
          file="K4K27_data_stat_maxK27over2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:77]

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
          file="K4K27_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27.over1, 
          file="K4K27_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27.over1, 
          file="K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)



##### Extraction of the genes with seasonality and overwriting the files #####

# Perform "09_edgeR.R", "10_Cosinor_model.R", and "12_Ampover1_edgeR_cosinor.R" in advance

K4 <- read.csv("K4K27_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv",header=T,sep=",")
K27 <- read.csv("K4K27_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv",header=T,sep=",")
K4.K27 <- merge(K4,K27,by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                            "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
K4.K27 <- K4.K27[,1:ncol(K4)]
colnames(K4.K27) <- gsub(".x", "", colnames(K4.K27))
colnames(K4.K27) <- gsub("mK", "maxK", colnames(K4.K27))

nrow(K4)
nrow(K27)
nrow(K4.K27)

write.csv(K4, 
          file="K4K27_data_stat_maxK4over2_ampK4over1.csv", 
          quote=F, row.names=F)
write.csv(K27, 
          file="K4K27_data_stat_maxK27over2_ampK27over1.csv", 
          quote=F, row.names=F)
write.csv(K4.K27, 
          file="K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv", 
          quote=F, row.names=F)



##### Extraction of genes with seasonality (K4 or K27 only) #####
K4.over1 <- read.csv("K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
K27.over1 <- read.csv("K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=T,sort=F)
K4.only.over1 <- K4.K27.over1[is.na(K4.K27.over1$K4_rep1_start.y),]
K27.only.over1 <- K4.K27.over1[is.na(K4.K27.over1$K4_rep1_start.x),]
K4.only.over1 <- K4.only.over1[,1:ncol(K4.over1)]
K27.only.over1 <- K27.only.over1[,c(1:7,(ncol(K4.over1)+1):ncol(K27.only.over1))]

colnames <- vector(length=ncol(K4.only.over1))
list <- strsplit(colnames(K4.only.over1), "\\.x")
for(i in 1:ncol(K4.only.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K4.only.over1) <- colnames

colnames <- vector(length=ncol(K27.only.over1))
list <- strsplit(colnames(K27.only.over1), "\\.y")
for(i in 1:ncol(K27.only.over1)){
  colnames[i] <- list[[i]][1]
}
colnames(K27.only.over1) <- colnames

sortlist <- order(K4.only.over1$ampK4, decreasing=T)
K4.only.over1 <- K4.only.over1[sortlist,]
sortlist <- order(K27.only.over1$ampK27, decreasing=T)
K27.only.over1 <- K27.only.over1[sortlist,]

nrow(K4.only.over1)
nrow(K27.only.over1)

write.csv(K4.only.over1, 
          file="K4K27_data_stat_maxK4over2_ampK4onlyover1.csv", 
          quote=F, row.names=F)
write.csv(K27.only.over1, 
          file="K4K27_data_stat_maxK27over2_ampK27onlyover1.csv", 
          quote=F, row.names=F)



##### Histogram of maximum #####
K4.K27 <- read.csv("K4K27_data_stat.csv",header=T,sep=",")

dir.create("../figs/histogram")

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


pdf("../figs/histogram/histogram_maxK4mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(maxK4.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,8), ylim=c(0,3000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,3000,1000),labels=F)
mtext(c("0","2","4","6","> 8"), side=1, at=seq(0,8,2),line=-0.3)
mtext(c("0","1,000","2,000","3,000"), side=2, at=seq(0,3000,1000),las=1,line=0.15)
mtext(expression(paste(log[2], "(maximum rpkm)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K4me3",side=3,line=0)
text(4.4,2500,"20,114",col="red")
arrows(2,0,2,3000,length=0,col="red")
arrows(2,2500,3,2500,length=0.05,col="red")
dev.off()

pdf("../figs/histogram/histogram_maxK27mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(maxK27.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,8), ylim=c(0,6000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,6000,2000),labels=F)
mtext(c("0","2","4","6","> 8"), side=1, at=seq(0,8,2),line=-0.3)
mtext(c("0","2,000","4,000","6,000"), side=2, at=seq(0,6000,2000),las=1,line=0.15)
mtext(expression(paste(log[2], "(maximum rpkm)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K27me3",side=3,line=0)
text(4.4,3000,"10,152",col="red")
arrows(2,0,2,6000,length=0,col="red")
arrows(2,3000,3,3000,length=0.05,col="red")
dev.off()



##### Histogram of amplitude #####
K4 <- read.csv("K4K27_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_data_stat_maxK27over2.csv",header=T,sep=",")
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

pdf("../figs/histogram/histogram_ampK4mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampK4.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,2), ylim=c(0,2000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,2000,1000),labels=F)
mtext(c("0","1","> 2"), side=1, at=seq(0,2,1),line=-0.3)
mtext(c("0","1,000","2,000"), side=2, at=seq(0,2000,1000),las=1,line=0.15)
mtext(expression(paste(log[2], "(seasonal amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K4me3",side=3,line=0)
text(1.5,1400,"2,346",col="red")
arrows(1,0,1,2000,length=0,col="red")
arrows(1,1400,1.2,1400,length=0.05,col="red")
dev.off()

pdf("../figs/histogram/histogram_ampK27mod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampK27.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,2), ylim=c(0,1000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,1000,500),labels=F)
mtext(c("0","1","> 2"), side=1, at=seq(0,2,1),line=-0.3)
mtext(c("0","500","1,000"), side=2, at=seq(0,1000,500),las=1,line=0.15)
mtext(expression(paste(log[2], "(seasonal amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("H3K27me3",side=3,line=0)
text(1.5,700,"1,911",col="red")
arrows(1,0,1,1000,length=0,col="red")
arrows(1,700,1.2,700,length=0.05,col="red")
dev.off()



# Venn diagram
dir.create("../figs/edgeR")
library(VennDiagram)
K4 <- read.csv("K4K27_data_stat_maxK4over2.csv",header=T,sep=",")
K27 <- read.csv("K4K27_data_stat_maxK27over2.csv",header=T,sep=",")

K4.over1 <- subset(K4,ampK4>1)
K27.over1 <- subset(K27,ampK27>1)
K4.K27.over1 <- 
  merge(K4.over1, K27.over1, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
K4.K27.over1 <- K4.K27.over1[,1:77]

nrow(K4.over1)
nrow(K27.over1)
nrow(K4.K27.over1)

pdf("../figs/edgeR/venn_K4K27_seasonalampK4K27over1.pdf",height=2.5,width=2.5)
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
K4 <- read.csv("K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
order <- order(K4$ampK4, decreasing=T)
K4.amp.order <- K4[order,]
write.csv(K4.amp.order, 
          file="K4K27_data_stat_maxK4over2_ampK4over1_ampK4order.csv",
          quote=F, row.names=F)



##### Ordering of genes according to ampK27 #####
K27 <- read.csv("K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
order <- order(K27$ampK27, decreasing=T)
K27.amp.order <- K27[order,]
write.csv(K27.amp.order, 
          file="K4K27_data_stat_maxK27over2_ampK27over1_ampK27order.csv",
          quote=F, row.names=F)



##### Calculation of monthly mean #####
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22],K4.K27.rep1[,11:22],K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37],K4.K27.rep1[,26:37],K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22],K4.K27.rep2[,11:22],K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37],K4.K27.rep2[,26:37],K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27.duplicate[,1:7]


K4.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K4.month) <- 
  c("K4.Nov", "K4.Dec", "K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", 
    "K4.May", "K4.Jun", "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct")
K27.month <- matrix(NA,nrow=nrow(K4),ncol=12)
colnames(K27.month) <- 
  c("K27.Nov", "K27.Dec", "K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", 
    "K27.May", "K27.Jun", "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct")

K4.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K4.min.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.max.date.fromJan1 <- rep(NA,length=nrow(K4))
K27.min.date.fromJan1 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(title, K4.month, K27.month, 
                K4.max.date.fromJan1, K4.min.date.fromJan1, 
                K27.max.date.fromJan1, K27.min.date.fromJan1)

for(i in 1:nrow(K4.K27)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  
  # monthly mean
  K4.K27$K4.Nov[i] <- mean(predK4$y[365<predK4$x & predK4$x<=365+30])
  K4.K27$K4.Dec[i] <- mean(predK4$y[365+30<predK4$x & predK4$x<=365+61])
  K4.K27$K4.Jan[i] <- mean(predK4$y[365+61<predK4$x & predK4$x<=365+92])
  K4.K27$K4.Feb[i] <- mean(predK4$y[365+92<predK4$x & predK4$x<=365+120])
  K4.K27$K4.Mar[i] <- mean(predK4$y[365+120<predK4$x & predK4$x<=365+151])
  K4.K27$K4.Apr[i] <- mean(predK4$y[365+151<predK4$x & predK4$x<=365+181])
  K4.K27$K4.May[i] <- mean(predK4$y[365+181<predK4$x & predK4$x<=365+212])
  K4.K27$K4.Jun[i] <- mean(predK4$y[365+212<predK4$x & predK4$x<=365+242])
  K4.K27$K4.Jul[i] <- mean(predK4$y[365+242<predK4$x & predK4$x<=365+273])
  K4.K27$K4.Aug[i] <- mean(predK4$y[365+273<predK4$x & predK4$x<=365+304])
  K4.K27$K4.Sep[i] <- mean(predK4$y[365+304<predK4$x & predK4$x<=365+334])
  K4.K27$K4.Oct[i] <- mean(predK4$y[365+334<predK4$x & predK4$x<=365+365])
  
  K4.K27$K27.Nov[i] <- mean(predK27$y[365<predK27$x & predK27$x<=365+30])
  K4.K27$K27.Dec[i] <- mean(predK27$y[365+30<predK27$x & predK27$x<=365+61])
  K4.K27$K27.Jan[i] <- mean(predK27$y[365+61<predK27$x & predK27$x<=365+92])
  K4.K27$K27.Feb[i] <- mean(predK27$y[365+92<predK27$x & predK27$x<=365+120])
  K4.K27$K27.Mar[i] <- mean(predK27$y[365+120<predK27$x & predK27$x<=365+151])
  K4.K27$K27.Apr[i] <- mean(predK27$y[365+151<predK27$x & predK27$x<=365+181])
  K4.K27$K27.May[i] <- mean(predK27$y[365+181<predK27$x & predK27$x<=365+212])
  K4.K27$K27.Jun[i] <- mean(predK27$y[365+212<predK27$x & predK27$x<=365+242])
  K4.K27$K27.Jul[i] <- mean(predK27$y[365+242<predK27$x & predK27$x<=365+273])
  K4.K27$K27.Aug[i] <- mean(predK27$y[365+273<predK27$x & predK27$x<=365+304])
  K4.K27$K27.Sep[i] <- mean(predK27$y[365+304<predK27$x & predK27$x<=365+334])
  K4.K27$K27.Oct[i] <- mean(predK27$y[365+334<predK27$x & predK27$x<=365+365])
  
  # Days from 1 Jan.
  if(max(predK4$y)>0){
    if((365+61 < predK4$x[predK4$y == max(predK4$y)]) &&
       (predK4$x[predK4$y == max(predK4$y)] <=365+365)){
      K4.K27$K4.max.date.fromJan1[i] <- 
        predK4$x[predK4$y == max(predK4$y)] - 365 -61
    }else{K4.K27$K4.max.date.fromJan1[i] <- 
      predK4$x[predK4$y == max(predK4$y)] + 365 - 365 -61}
    if((365+61 < predK4$x[predK4$y == min(predK4$y)]) && 
       (predK4$x[predK4$y == min(predK4$y)] <=365+365)){
      K4.K27$K4.min.date.fromJan1[i] <- 
        predK4$x[predK4$y == min(predK4$y)] - 365 -61
    }else{K4.K27$K4.min.date.fromJan1[i] <- 
      predK4$x[predK4$y == min(predK4$y)] + 365 - 365 -61}	
  }else{}
  
  if(max(predK27$y)>0){
    if((365+61 < predK27$x[predK27$y == max(predK27$y)]) && 
       (predK27$x[predK27$y == max(predK27$y)] <=365+365)){
      K4.K27$K27.max.date.fromJan1[i] <- 
        predK27$x[predK27$y == max(predK27$y)] - 365 -61
    }else{K4.K27$K27.max.date.fromJan1[i] <- 
      predK27$x[predK27$y == max(predK27$y)] + 365 - 365 -61}
    if((365+61 < predK27$x[predK27$y == min(predK27$y)]) && 
       (predK27$x[predK27$y == min(predK27$y)] <=365+365)){
      K4.K27$K27.min.date.fromJan1[i] <- 
        predK27$x[predK27$y == min(predK27$y)] - 365 -61
    }else{K4.K27$K27.min.date.fromJan1[i] <- 
      predK27$x[predK27$y == min(predK27$y)] + 365 - 365 -61}
  }else{}
}

K4.K27 <- read.csv("K4K27_monthmean.csv",header=T,sep=",")
K4.K27.2 <- K4.K27[,c(1:7,10:19,8,9,22:31,20,21,32:35)]

write.csv(K4.K27.2, file="K4K27_monthmean.csv", quote=F, row.names=F)



#####  Extraction of genes with amp > 1 #####
K4 <- read.csv("K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
K27 <- read.csv("K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
K4K27<-read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv",header=T,sep=",")
K4only<-read.csv("K4K27_data_stat_maxK4over2_ampK4onlyover1.csv",header=T,sep=",")
K27only<-read.csv("K4K27_data_stat_maxK27over2_ampK27onlyover1.csv",header=T,sep=",")

K4K27mean <- read.csv("K4K27_monthmean.csv",header=T,sep=",")

K4.seasonality.mean <- 
  merge(K4K27mean,K4, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4.seasonality.mean <- K4.seasonality.mean[,1:ncol(K4K27mean)]

K27.seasonality.mean <- 
  merge(K4K27mean,K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27.seasonality.mean <- K27.seasonality.mean[,1:ncol(K4K27mean)]

K4K27.seasonality.mean <- 
  merge(K4K27mean,K4K27, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4K27.seasonality.mean <- K4K27.seasonality.mean[,1:ncol(K4K27mean)]

K4only.seasonality.mean <- 
  merge(K4K27mean,K4only, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4only.seasonality.mean <- K4only.seasonality.mean[,1:ncol(K4K27mean)]

K27only.seasonality.mean <- 
  merge(K4K27mean,K27only, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K27only.seasonality.mean <- K27only.seasonality.mean[,1:ncol(K4K27mean)]

nrow(K4.seasonality.mean)
nrow(K27.seasonality.mean)
nrow(K4K27.seasonality.mean)
nrow(K4only.seasonality.mean)
nrow(K27only.seasonality.mean)

write.csv(K4.seasonality.mean, 
          file="K4K27_monthmean_ampK4over1.csv", quote=F, row.names=F)
write.csv(K27.seasonality.mean, 
          file="K4K27_monthmean_ampK27over1.csv", quote=F, row.names=F)
write.csv(K4K27.seasonality.mean, 
          file="K4K27_monthmean_ampK4K27over1.csv", quote=F, row.names=F)
write.csv(K4only.seasonality.mean, 
          file="K4K27_monthmean_ampK4onlyover1.csv", quote=F, row.names=F)
write.csv(K27only.seasonality.mean, 
          file="K4K27_monthmean_ampK27onlyover1.csv", quote=F, row.names=F)



##### Correlation between monthly means　(amp over1) #############################
K4 <- read.csv("K4K27_monthmean_ampK4over1.csv",header=T,sep=",")
K4.2 <- K4[,c(9:19,8)]
K27 <- read.csv("K4K27_monthmean_ampK27over1.csv",header=T,sep=",")
K27.2 <- K27[,c(21:31,20)]

names <- c("2","4","6","8","10","12")

dir.create("../figs/heatmap")
pdf("../figs/heatmap/CorMonth_K4_monthmean_ampK4over1.pdf",width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(mgp=c(0,0.3,0))
par(ps=6)
par(xpd=T)
image(cor(K4.2,method="spearman"),
	  zlim=c(min(c(cor(K4.2,method="spearman"),cor(K27.2,method="spearman"))),1),
	  col=ng.po.colors(64),axes=F)
axis(side=1,at=seq(0,1-1/12,length=6),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,1-1/12,length=6),las=1,tcl=-0.1,labels=F)
mtext(names,at=seq(0,1-1/12,length=6),side=1,line=-0.32,las=1)
mtext(names,at=seq(0,1-1/12,length=6),side=2,line=0.15,las=1)
box()
mtext("Month", side=1, line=0.1)
mtext("Month", side=2, line=0.5)
mtext("H3K4me3", side=3, line=0)
image.plot(legend.only=T,
		   zlim=c(min(c(cor(K4.2,method="spearman"),cor(K27.2,method="spearman"))),1),
		   legend.width=0.3,tcl=-0.1,col=ng.po.colors(64),legend.mar=3.4)
mtext("Correlation",side=4,line=1.7)
dev.off()


pdf("../figs/heatmap/CorMonth_K27_monthmean_ampK27over1.pdf",width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(mgp=c(0,0.3,0))
par(ps=6)
par(xpd=T)
image(cor(K27.2,method="spearman"),
	  zlim=c(min(c(cor(K4.2,method="spearman"),cor(K27.2,method="spearman"))),1),
	  col=ng.po.colors(64),axes=F)
axis(side=1,at=seq(0,1-1/12,length=6),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,1-1/12,length=6),las=1,tcl=-0.1,labels=F)
mtext(names,at=seq(0,1-1/12,length=6),side=1,line=-0.32,las=1)
mtext(names,at=seq(0,1-1/12,length=6),side=2,line=0.15,las=1)
box()
mtext("Month", side=1, line=0.1)
mtext("Month", side=2, line=0.5)
mtext("H3K27me3", side=3, line=0)
image.plot(legend.only=T,
		   zlim=c(min(c(cor(K4.2,method="spearman"),cor(K27.2,method="spearman"))),1),
		   legend.width=0.3,tcl=-0.1,col=ng.po.colors(64),legend.mar=3.4)
mtext("Correlation",side=4,line=1.7)
dev.off()



# Line graph  Feb, May, Aug, Nov and each month
K4 <- read.csv("K4K27_monthmean_ampK4over1.csv",header=T,sep=",")
K4.2 <- K4[,c(9:19,8)]
K27 <- read.csv("K4K27_monthmean_ampK27over1.csv",header=T,sep=",")
K27.2 <- K27[,c(21:31,20)]

names <- c("2","4","6","8","10","12")

# K4 
cor.vectK4.Feb <- vector(length=12)
conf.low.vectK4.Feb <- vector(length=12)
conf.high.vectK4.Feb <- vector(length=12)
cor.vectK4.May <- vector(length=12)
conf.low.vectK4.May <- vector(length=12)
conf.high.vectK4.May <- vector(length=12)
cor.vectK4.Aug <- vector(length=12)
conf.low.vectK4.Aug <- vector(length=12)
conf.high.vectK4.Aug <- vector(length=12)
cor.vectK4.Nov <- vector(length=12)
conf.low.vectK4.Nov <- vector(length=12)
conf.high.vectK4.Nov <- vector(length=12)

for(i in 1:12){
	cor.vectK4.Feb[i] <- cor.test(rank(K4.2[,1]),rank(K4.2[,i]))[[4]]
	conf.low.vectK4.Feb[i] <- cor.test(rank(K4.2[,1]),rank(K4.2[,i]))[[9]][1]
	conf.high.vectK4.Feb[i] <- cor.test(rank(K4.2[,1]),rank(K4.2[,i]))[[9]][2]
	cor.vectK4.May[i] <- cor.test(rank(K4.2[,4]),rank(K4.2[,i]))[[4]]
	conf.low.vectK4.May[i] <- cor.test(rank(K4.2[,4]),rank(K4.2[,i]))[[9]][1]
	conf.high.vectK4.May[i] <- cor.test(rank(K4.2[,4]),rank(K4.2[,i]))[[9]][2]
	cor.vectK4.Aug[i] <- cor.test(rank(K4.2[,7]),rank(K4.2[,i]))[[4]]
	conf.low.vectK4.Aug[i] <- cor.test(rank(K4.2[,7]),rank(K4.2[,i]))[[9]][1]
	conf.high.vectK4.Aug[i] <- cor.test(rank(K4.2[,7]),rank(K4.2[,i]))[[9]][2]
	cor.vectK4.Nov[i] <- cor.test(rank(K4.2[,10]),rank(K4.2[,i]))[[4]]
	conf.low.vectK4.Nov[i] <- cor.test(rank(K4.2[,10]),rank(K4.2[,i]))[[9]][1]
	conf.high.vectK4.Nov[i] <- cor.test(rank(K4.2[,10]),rank(K4.2[,i]))[[9]][2]
}

# K27 
cor.vectK27.Feb <- vector(length=12)
conf.low.vectK27.Feb <- vector(length=12)
conf.high.vectK27.Feb <- vector(length=12)
cor.vectK27.May <- vector(length=12)
conf.low.vectK27.May <- vector(length=12)
conf.high.vectK27.May <- vector(length=12)
cor.vectK27.Aug <- vector(length=12)
conf.low.vectK27.Aug <- vector(length=12)
conf.high.vectK27.Aug <- vector(length=12)
cor.vectK27.Nov <- vector(length=12)
conf.low.vectK27.Nov <- vector(length=12)
conf.high.vectK27.Nov <- vector(length=12)

for(i in 1:12){
	cor.vectK27.Feb[i] <- cor.test(rank(K27.2[,1]),rank(K27.2[,i]))[[4]]
	conf.low.vectK27.Feb[i] <- cor.test(rank(K27.2[,1]),rank(K27.2[,i]))[[9]][1]
	conf.high.vectK27.Feb[i] <- cor.test(rank(K27.2[,1]),rank(K27.2[,i]))[[9]][2]
	cor.vectK27.May[i] <- cor.test(rank(K27.2[,4]),rank(K27.2[,i]))[[4]]
	conf.low.vectK27.May[i] <- cor.test(rank(K27.2[,4]),rank(K27.2[,i]))[[9]][1]
	conf.high.vectK27.May[i] <- cor.test(rank(K27.2[,4]),rank(K27.2[,i]))[[9]][2]
	cor.vectK27.Aug[i] <- cor.test(rank(K27.2[,7]),rank(K27.2[,i]))[[4]]
	conf.low.vectK27.Aug[i] <- cor.test(rank(K27.2[,7]),rank(K27.2[,i]))[[9]][1]
	conf.high.vectK27.Aug[i] <- cor.test(rank(K27.2[,7]),rank(K27.2[,i]))[[9]][2]
	cor.vectK27.Nov[i] <- cor.test(rank(K27.2[,10]),rank(K27.2[,i]))[[4]]
	conf.low.vectK27.Nov[i] <- cor.test(rank(K27.2[,10]),rank(K27.2[,i]))[[9]][1]
	conf.high.vectK27.Nov[i] <- cor.test(rank(K27.2[,10]),rank(K27.2[,i]))[[9]][2]
}


pdf("../figs/heatmap/CorFebMayAugNov_K4K27_monthmean.pdf", height=1.1,width=2)
  set.panel(2,2,relax=T)
  par(oma=c(0.4,0.3,0,0))
  par(mar=c(0.6,1,0.2,0.1))
  par(mgp=c(0,0.15,0))
  par(ps=6)
  par(cex=1)
  par(xpd=T)
  
# February
  plot(1:12,cor.vectK4.Feb,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.5,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  arrows(1:12, conf.low.vectK4.Feb, 1:12, conf.high.vectK4.Feb, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0.9,0,0))

  par(new=T)
  plot(1:12,cor.vectK27.Feb,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.5,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  arrows(1:12, conf.low.vectK27.Feb, 1:12, conf.high.vectK27.Feb, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0,0,0.7))

  axis(side=1,at=seq(1,11,2),label=F,tcl=-0.1)
  mtext(names,side=1,at=seq(1,11,2),line=-0.35)
  #mtext("Month",side=1,line=0.05)
  axis(side=2,at=seq(0.6,1,0.2),tcl=-0.1,las=1)
  mtext("Correlation",side=2,line=0.65)
  text(6.5,0.96,"February")
  #mtext("February and each month",side=3,line=0)
  #text(c(4,9,10),conf.high.vect[c(4,9,10)]+0.025,"*")
  box()
  

# May
  plot(1:12,cor.vectK4.May,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.65,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  arrows(1:12, conf.low.vectK4.May, 1:12, conf.high.vectK4.May, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0.9,0,0))

  par(new=T)
  plot(1:12,cor.vectK27.May,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.65,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  arrows(1:12, conf.low.vectK27.May, 1:12, conf.high.vectK27.May, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0,0,0.7))

  axis(side=1,at=seq(1,11,2),label=F,tcl=-0.1)
  mtext(names,side=1,at=seq(1,11,2),line=-0.35)
  #mtext("Month",side=1,line=0.05)
  axis(side=2,at=seq(0.7,1,0.1),tcl=-0.1,las=1)
  #mtext("Correlation",side=2,line=0.65)
  text(6.5,0.71,"May")
  #mtext("May and each month",side=3,line=0)
  #text(c(4,9,10),conf.high.vect[c(4,9,10)]+0.025,"*")
  box()


# August
  plot(1:12,cor.vectK4.Aug,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.5,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  arrows(1:12, conf.low.vectK4.Aug, 1:12, conf.high.vectK4.Aug, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0.9,0,0))

  par(new=T)
  plot(1:12,cor.vectK27.Aug,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.5,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  arrows(1:12, conf.low.vectK27.Aug, 1:12, conf.high.vectK27.Aug, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0,0,0.7))

  axis(side=1,at=seq(1,11,2),label=F,tcl=-0.1)
  mtext(names,side=1,at=seq(1,11,2),line=-0.35)
  mtext("Month",side=1,line=0.05)
  axis(side=2,at=seq(0.6,1,0.2),tcl=-0.1,las=1)
  mtext("Correlation",side=2,line=0.65)
  text(6.5,0.58,"August")
  #mtext("August and each month",side=3,line=0)
  #text(c(4,9,10),conf.high.vect[c(4,9,10)]+0.025,"*")
  box()


# November
  plot(1:12,cor.vectK4.Nov,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.75,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  arrows(1:12, conf.low.vectK4.Nov, 1:12, conf.high.vectK4.Nov, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0.9,0,0))

  par(new=T)
  plot(1:12,cor.vectK27.Nov,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(0.75,1.05),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  arrows(1:12, conf.low.vectK27.Nov, 1:12, conf.high.vectK27.Nov, 
  		 length=.01, angle=90, code=3, lwd=0.5, col=rgb(0,0,0.7))

  axis(side=1,at=seq(1,11,2),label=F,tcl=-0.1)
  mtext(names,side=1,at=seq(1,11,2),line=-0.35)
  mtext("Month",side=1,line=0.05)
  axis(side=2,at=seq(0.8,1,0.1),tcl=-0.1,las=1)
  #mtext("Correlation",side=2,line=0.65)
  text(5,1.01,"November")
  #mtext("November and each month",side=3,line=0)
  #text(c(4,9,10),conf.high.vect[c(4,9,10)]+0.025,"*")
  box()
  
dev.off()	



pdf("../figs/heatmap/CorFebMayAugNov_K4K27_monthmean_legend.pdf",width=1.2,height=0.15)
	par(mar=c(0,0,0,0))
	par(xpd=T)
	par(ps=6)
	plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F, ann=F)
	legend(-0.05,1.6,legend="H3K27me3",col=rgb(0,0,0.7),
		   lty=1,pch=16,pt.cex=0.3,bty="n",x.intersp=0.2,y.intersp=0.5,seg.len=0.7)
	legend(0.53,1.6,legend="H3K4me3",col=rgb(0.9,0,0),
		   lty=1,pch=16,pt.cex=0.3,bty="n",x.intersp=0.2,y.intersp=0.5,seg.len=0.7)
dev.off()



##### Correlation between monthly means (between K4 and K27) #####
K4K27 <- read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")
K4 <- K4K27[,c(9:19,8)]
K27 <- K4K27[,c(21:31,20)]

cor.mat <- matrix(NA,nrow=12,ncol=12)
for(i in 1:12){
	for(j in 1:12){
  		cor.mat[i,j] <- cor(K4[,i],K27[,j],method="spearman")
  	}
}

names <- c("2","4","6","8","10","12")

pdf("../figs/heatmap/CorMonth_K4K27_monthmean_ampK4K27over1.pdf",width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(mgp=c(0,0.3,0))
par(ps=6)
par(xpd=T)
image(cor.mat,zlim=c(min(cor.mat),max(cor.mat)),col=ng.po.colors(64),axes=F)
axis(side=1,at=seq(0,1-1/12,length=6),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,1-1/12,length=6),las=1,tcl=-0.1,labels=F)
mtext(names,at=seq(0,1-1/12,length=6),side=1,line=-0.32,las=1)
mtext(names,at=seq(0,1-1/12,length=6),side=2,line=0.15,las=1)
box()
mtext("Month (H3K4me3)", side=1, line=0.1)
mtext("Month (H3K27me3)", side=2, line=0.5)
#mtext("H3K4me3", side=3, line=0)
image.plot(legend.only=T,zlim=c(min(cor.mat),max(cor.mat)),legend.width=0.3,tcl=-0.1,
			col=ng.po.colors(64),legend.mar=3.4)
mtext("Correlation",side=4,line=1.9)
dev.off()



##### Correlation at the same month (between K4 and K27) #####
K4K27 <- read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")
K4 <- K4K27[,c(8:19)]
K27 <- K4K27[,c(20:31)]
names <- c("1","3","5","7","9","11")

cor.vect <- vector(length=12)
for(i in 1:12){
  		cor.vect[i] <- cor.test(rank(K4[,i]),rank(K27[,i]))[[4]]
}

pdf("../figs/heatmap/CorSameMonth_K4K27_monthmean_ampK4K27over1.pdf", height=1,width=1.4)
  par(oma=c(0,0,0,0))
  par(mar=c(1,2,1.4,0.6))
  par(mgp=c(0,0.3,0))
  par(ps=6)
  
  plot(1:12,cor.vect,
       type="o",pch=16,las=1,tcl=-0.2,ylim=c(-0.47,-0.25),
       xlab="",ylab="",axes=F, ann=F, cex=0.3)

  axis(side=1,at=seq(1,11,2),label=F,tcl=-0.1)
  mtext(names,side=1,at=seq(1,11,2),line=-0.35)
  mtext("Month",side=1,line=0.05)
  axis(side=2,at=seq(-0.4,-0.3,0.1),label=F,tcl=-0.1)
  mtext(c("-0.4","-0.3"),side=2,
  		at=seq(-0.4,-0.3,0.1),line=0.15,las=1)
  #mtext("Correlation",side=2,line=1.1)
  mtext("Correlation between",side=3,line=0.8)
  mtext("monthly means of",side=3,line=0.4)
  mtext("H3K27me3 and H3K4me3",side=3,line=0)
  text(1:2,cor.vect[1:2]+0.025,"*")
  box()
dev.off()	



##### Statictical test of cor differences #####
for(i in 1){
	Nov <- as.data.frame(cbind(rank(K4[,11]),rank(K27[,11])))
	colnames(Nov) <- c("K4.Nov","K27.Nov")
	Mon <- as.data.frame(cbind(rank(K4[,i]),rank(K27[,i])))
	colnames(Mon) <- c("K4.Mon","K27.Mon")
	Nov.Mon.list <- list(Nov,Mon)
	names(Nov.Mon.list) <- c("Nov","Mon")
}

cocor(~K4.Nov + K27.Nov | K4.Mon + K27.Mon, Nov.Mon.list)
#significant: 1 (Jan), 2 (Feb)



##### Spearman's rank correlation between seasonal changes of K27 and K4 #####
K4K27<-read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")
K4 <- K4K27[,8:19]
K27 <- K4K27[,20:31]
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
title <- K4K27[,1:7]

Cor.K4.K27 <- rep(NA,length=nrow(K4))
cor.mat <- cbind(title, Cor.K4.K27)
for(i in 1:nrow(K4)){
  cor.mat[i,8] <- round(cor(K4[i,],K27[i,], method="spearman"), digits=2)
}

sortlist <- order(as.numeric(cor.mat[,8]))
cor.mat.sorted <- cor.mat[sortlist,]
write.csv(cor.mat.sorted, 
          file="K4K27_ampK4K27over1_cor.csv",quote=F, row.names=F)

pdf("../figs/histogram/histogram_K4K27_ampK4K27over1_cor.pdf",width=1.6,height=1.2)
par(ps=6)
par(mar=c(2,2,1.3,1))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(as.numeric(cor.mat[,8]), las=1, tcl=-0.2, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-1,1), ylim=c(0,150),breaks=24,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0,0),labels=F)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,150,50))
label <- c("-1.0","-0.5","0.0","0.5","1.0")
mtext(label,side=1,at=c(-1,-0.5,0,0.5,1),line=-0.35,las=1)
mtext("Correlation between",side=1,line=0.1)
mtext("seasonal dynamics of",side=1,line=0.5)
mtext("H3K27me3 and H3K4me3",side=1,line=0.9)
mtext("No. of genes",side=2,line=0.8)
#mtext("Genes with",side=3,line=0.4)
#mtext("dual-seasonality",side=3,line=0)
arrows(0,0,0,150,length=0,col="black",lty="dashed",lwd=0.5)
dev.off()



##### Histogram of Peak month #####

# Determination of peak month
K4 <- read.csv("K4K27_monthmean_ampK4over1.csv",header=T,sep=",")
K27 <- read.csv("K4K27_monthmean_ampK27over1.csv",header=T,sep=",")

K4.maxK4 <- rep(NA, length=nrow(K4))
K4.minK4 <- rep(NA, length=nrow(K4))
K4.maxK27 <- rep(NA, length=nrow(K4))
K4.minK27 <- rep(NA, length=nrow(K4))
K27.maxK4 <- rep(NA, length=nrow(K27))
K27.minK4 <- rep(NA, length=nrow(K27))
K27.maxK27 <- rep(NA, length=nrow(K27))
K27.minK27 <- rep(NA, length=nrow(K27))

for(i in 1:nrow(K4)){
  K4.maxK4[i] <- as.numeric(which.max(K4[i,8:19]))
  K4.minK4[i] <- as.numeric(which.min(K4[i,8:19]))
  K4.maxK27[i] <- as.numeric(which.max(K4[i,20:31]))
  K4.minK27[i] <- as.numeric(which.min(K4[i,20:31]))
}

for(i in 1:nrow(K27)){
  K27.maxK4[i] <- as.numeric(which.max(K27[i,8:19]))
  K27.minK4[i] <- as.numeric(which.min(K27[i,8:19]))
  K27.maxK27[i] <- as.numeric(which.max(K27[i,20:31]))
  K27.minK27[i] <- as.numeric(which.min(K27[i,20:31]))
}

K4 <- cbind(K4, K4.maxK4, K4.minK4, K4.maxK27, K4.minK27)
K27 <- cbind(K27, K27.maxK4, K27.minK4, K27.maxK27, K27.minK27)

write.csv(K4, file="K4K27_monthmean_ampK4over1_peakmonth.csv", 
          quote=F, row.names=F)
write.csv(K27, file="K4K27_monthmean_ampK27over1_peakmonth.csv", 
          quote=F, row.names=F)


# Count of peak month
K4 <- read.csv("K4K27_monthmean_ampK4over1_peakmonth.csv",header=T,sep=",")
K27 <- read.csv("K4K27_monthmean_ampK27over1_peakmonth.csv",header=T,sep=",")
maxK4<-K4$K4.maxK4
minK4<-K4$K4.minK4
maxK27<-K27$K27.maxK27
minK27<-K27$K27.minK27
maxK4count <- table(maxK4)
minK4count <- table(minK4)
maxK27count <- table(maxK27)
minK27count <- table(minK27)
month <- 1:12
count <- cbind(month, maxK4count, minK4count, maxK27count, minK27count)
write.csv(count, file="peakmonth_count.csv", quote=F, row.names=F)

# Histogram
z <- read.csv("peakmonth_count.csv",header=T,sep=",")
label <- c("1","","3","","5","","7","","9","","11","")

pdf("../figs/histogram/histogram_maxK4.pdf", width=1.14, height=0.75)
par(mar=c(1.5,2.5,1.2,0.2))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)
barplot(z$maxK4count, beside=T,ylim=c(0,500), las=1, tcl=-0.1, 
        axes=F,ann=F, mgp=c(3,0,0.5),space=0)
#axis(side=1, at=seq(0,12,1),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,500,250),tcl=-0.1,las=2,pos=-0.5)
mtext(label,side=1,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),line=-0.35)
mtext("Month", side=1, at=-3.5, line=-0.35)
mtext("No. of", side=2, line=1.15)
mtext("genes", side=2, line=0.8)
mtext("Peak month", side=3, line=0)
mtext("H3K4me3", at=-3.6, col=rgb(0.9,0,0), side=3, line=0.4)
dev.off()

pdf("../figs/histogram/histogram_minK4.pdf", width=1.14, height=0.75)
par(mar=c(1.5,2.5,1.2,0.2))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)
barplot(z$minK4count, beside=T,ylim=c(0,1000), las=1, tcl=-0.1, 
        axes=F,ann=F, mgp=c(3,0,0.5),space=0)
#axis(side=1, at=seq(0,12,1),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,1000,500),tcl=-0.1,las=2,pos=-0.5)
mtext(label,side=1,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),line=-0.35)
mtext("Month", side=1, at=-3.5, line=-0.35)
mtext("No. of", side=2, line=1.25)
mtext("genes", side=2, line=0.9)
mtext("Bottom month", side=3, line=0)
mtext("H3K4me3", at=-3.5, col=rgb(0.9,0,0), side=3, line=0.4)
dev.off()

pdf("../figs/histogram/histogram_maxK27.pdf", width=1.14, height=0.75)
par(mar=c(1.5,2.5,1.2,0.2))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)
barplot(z$maxK27count, beside=T,ylim=c(0,500), las=1, tcl=-0.1, 
        axes=F,ann=F, mgp=c(3,0,0.5),space=0)
#axis(side=1, at=seq(0,12,1),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,500,250),tcl=-0.1,las=2,pos=-0.5)
mtext(label,side=1,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),line=-0.35)
mtext("Month", side=1, at=-3.5, line=-0.35)
mtext("No. of", side=2, line=1.15)
mtext("genes", side=2, line=0.8)
mtext("Peak month", side=3, line=0)
mtext("H3K27me3", at=-3.1, col=rgb(0,0,0.7), side=3, line=0.4)
dev.off()

pdf("../figs/histogram/histogram_minK27.pdf", width=1.14, height=0.75)
par(mar=c(1.5,2.5,1.2,0.2))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)
barplot(z$minK27count, beside=T,ylim=c(0,1000), las=1, tcl=-0.1, 
        axes=F,ann=F, mgp=c(3,0,0.5),space=0)
#axis(side=1, at=seq(0,12,1),tcl=-0.1,labels=F)
axis(side=2,at=seq(0,1000,500),tcl=-0.1,las=2,pos=-0.5)
mtext(label,side=1,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),line=-0.35)
mtext("Month", side=1, at=-3.5, line=-0.35)
mtext("No. of", side=2, line=1.25)
mtext("genes", side=2, line=0.9)
mtext("Bottom month", side=3, line=0)
mtext("H3K27me3", at=-3.5, col=rgb(0,0,0.7), side=3, line=0.4)
dev.off()



##### Heatmap of monthly mean #####
# sorting of genes according to the date of maxK4 (ampK4over1)
y<-read.csv("K4K27_monthmean_ampK4over1.csv",header=T,sep=",")

sortlist <- order(y$K4.max.date.fromJan1)
y2 <- y[sortlist,]

write.csv(y2, file="K4K27_monthmean_ampK4over1_maxK4dateorder.csv", 
          quote=F, row.names=F)

# Normalization
K4<-read.csv("K4K27_monthmean_ampK4over1_maxK4dateorder.csv",
             header=T,sep=",")

K4.2 <- K4[,8:19]
K4.2 <- as.matrix(K4.2)
K4.normalized <- matrix(NA,nrow=nrow(K4.2),ncol=12)
for(i in 1:nrow(K4.2)){
  for(j in 1:12){
    K4.normalized[i,j] <- 
      (K4.2[i,j]-min(K4.2[i,]))/(max(K4.2[i,])-min(K4.2[i,]))
  }
}
colnames(K4.normalized) <- 
  c("K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", "K4.May", "K4.Jun",
    "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct", "K4.Nov", "K4.Dec")

K27.2 <- K4[,20:31]
K27.2 <- as.matrix(K27.2)
K27.normalized <- matrix(NA,nrow=nrow(K27.2),ncol=12)
for(i in 1:nrow(K27.2)){
  for(j in 1:12){
    K27.normalized[i,j] <- 
      (K27.2[i,j]-min(K27.2[i,]))/(max(K27.2[i,])-min(K27.2[i,]))
  }
}
colnames(K27.normalized) <- 
  c("K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", "K27.May", "K27.Jun", 
    "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct", "K27.Nov", "K27.Dec")

K4K27.normalized <- cbind(K4[,1:7],K4.normalized, K27.normalized, K4[,32:33])

write.csv(K4K27.normalized, 
          file="K4K27_monthmean_ampK4over1_maxK4dateorder_normalized.csv", 
          quote=F, row.names=F)

# Heatmap
dir.create("../figs/heatmap")

K4<-read.csv("K4K27_monthmean_ampK4over1_maxK4dateorder_normalized.csv",
             header=T,sep=",")
K27 <- t(apply(K4[,20:31],2,rev))
K4 <- t(apply(K4[,8:19],2,rev))

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK4over1_maxK4dateorder_normalized.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0.9,0,0)))
image.plot(K4,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K4me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK4over1_maxK4dateorder_normalized_K27.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0,0,0.7)))
image.plot(K27,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K27me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()


# sorting of genes according to the date of maxK27 (ampK27over1)
y<-read.csv("K4K27_monthmean_ampK27over1.csv",header=T,sep=",")

sortlist <- order(y$K27.max.date.fromJan1)
y2 <- y[sortlist,]

write.csv(y2, file="K4K27_monthmean_ampK27over1_maxK27dateorder.csv", 
          quote=F, row.names=F)

# Normalization
K27<-read.csv("K4K27_monthmean_ampK27over1_maxK27dateorder.csv",
              header=T,sep=",")

K4.2 <- K27[,8:19]
K4.2 <- as.matrix(K4.2)
K4.normalized <- matrix(NA,nrow=nrow(K4.2),ncol=12)
for(i in 1:nrow(K4.2)){
  for(j in 1:12){
    K4.normalized[i,j] <- 
      (K4.2[i,j]-min(K4.2[i,]))/(max(K4.2[i,])-min(K4.2[i,]))
  }
}
colnames(K4.normalized) <- 
  c("K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", "K4.May", "K4.Jun", 
    "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct", "K4.Nov", "K4.Dec")

K27.2 <- K27[,20:31]
K27.2 <- as.matrix(K27.2)
K27.normalized <- matrix(NA,nrow=nrow(K27.2),ncol=12)
for(i in 1:nrow(K27.2)){
  for(j in 1:12){
    K27.normalized[i,j] <- 
      (K27.2[i,j]-min(K27.2[i,]))/(max(K27.2[i,])-min(K27.2[i,]))
  }
}
colnames(K27.normalized) <- 
  c("K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", "K27.May", "K27.Jun", 
    "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct", "K27.Nov", "K27.Dec")

K4K27.normalized <- cbind(K27[,1:7],K4.normalized, K27.normalized, K27[,34:35])

write.csv(K4K27.normalized, 
          file="K4K27_monthmean_ampK27over1_maxK27dateorder_normalized.csv", 
          quote=F, row.names=F)

# Heatmap
K27<-read.csv("K4K27_monthmean_ampK27over1_maxK27dateorder_normalized.csv",
              header=T,sep=",")
K4 <- t(apply(K27[,8:19],2,rev))
K27 <- t(apply(K27[,20:31],2,rev))

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK27over1_maxK27dateorder_normalized.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0,0,0.7)))
image.plot(K27,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K27me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK27over1_maxK27dateorder_normalized_K4.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0.9,0,0)))
image.plot(K4,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K4me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()



##### Heatmap of monthly mean for dual seasonality genes #####

# sorting of genes according to the date of maxK4 (ampK4over1)
y<-read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")

sortlist <- order(y$K4.max.date.fromJan1)
y2 <- y[sortlist,]

write.csv(y2, file="K4K27_monthmean_ampK4K27over1_maxK4dateorder.csv", 
          quote=F, row.names=F)

# Normalization
K4<-read.csv("K4K27_monthmean_ampK4K27over1_maxK4dateorder.csv",
             header=T,sep=",")

K4.2 <- K4[,8:19]
K4.2 <- as.matrix(K4.2)
K4.normalized <- matrix(NA,nrow=nrow(K4.2),ncol=12)
for(i in 1:nrow(K4.2)){
  for(j in 1:12){
    K4.normalized[i,j] <- 
      (K4.2[i,j]-min(K4.2[i,]))/(max(K4.2[i,])-min(K4.2[i,]))
  }
}
colnames(K4.normalized) <- 
  c("K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", "K4.May", "K4.Jun",
    "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct", "K4.Nov", "K4.Dec")

K27.2 <- K4[,20:31]
K27.2 <- as.matrix(K27.2)
K27.normalized <- matrix(NA,nrow=nrow(K27.2),ncol=12)
for(i in 1:nrow(K27.2)){
  for(j in 1:12){
    K27.normalized[i,j] <- 
      (K27.2[i,j]-min(K27.2[i,]))/(max(K27.2[i,])-min(K27.2[i,]))
  }
}
colnames(K27.normalized) <- 
  c("K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", "K27.May", "K27.Jun", 
    "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct", "K27.Nov", "K27.Dec")

K4K27.normalized <- cbind(K4[,1:7],K4.normalized, K27.normalized, K4[,32:33])

write.csv(K4K27.normalized, 
          file="K4K27_monthmean_ampK4K27over1_maxK4dateorder_normalized.csv", 
          quote=F, row.names=F)

# Heatmap
K4<-read.csv("K4K27_monthmean_ampK4K27over1_maxK4dateorder_normalized.csv",
             header=T,sep=",")
K27 <- t(apply(K4[,20:31],2,rev))
K4 <- t(apply(K4[,8:19],2,rev))

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK4K27over1_maxK4dateorder_normalized.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0.9,0,0)))
image.plot(K4,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K4me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK4K27over1_maxK4dateorder_normalized_K27.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0,0,0.7)))
image.plot(K27,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K27me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()


# sorting of genes according to the date of maxK27 (ampK27over1)
y<-read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")

sortlist <- order(y$K27.max.date.fromJan1)
y2 <- y[sortlist,]

write.csv(y2, file="K4K27_monthmean_ampK4K27over1_maxK27dateorder.csv", 
          quote=F, row.names=F)

# Normalization
K27<-read.csv("K4K27_monthmean_ampK4K27over1_maxK27dateorder.csv",
              header=T,sep=",")

K4.2 <- K27[,8:19]
K4.2 <- as.matrix(K4.2)
K4.normalized <- matrix(NA,nrow=nrow(K4.2),ncol=12)
for(i in 1:nrow(K4.2)){
  for(j in 1:12){
    K4.normalized[i,j] <- 
      (K4.2[i,j]-min(K4.2[i,]))/(max(K4.2[i,])-min(K4.2[i,]))
  }
}
colnames(K4.normalized) <- 
  c("K4.Jan", "K4.Feb", "K4.Mar", "K4.Apr", "K4.May", "K4.Jun", 
    "K4.Jul", "K4.Aug", "K4.Sep", "K4.Oct", "K4.Nov", "K4.Dec")

K27.2 <- K27[,20:31]
K27.2 <- as.matrix(K27.2)
K27.normalized <- matrix(NA,nrow=nrow(K27.2),ncol=12)
for(i in 1:nrow(K27.2)){
  for(j in 1:12){
    K27.normalized[i,j] <- 
      (K27.2[i,j]-min(K27.2[i,]))/(max(K27.2[i,])-min(K27.2[i,]))
  }
}
colnames(K27.normalized) <- 
  c("K27.Jan", "K27.Feb", "K27.Mar", "K27.Apr", "K27.May", "K27.Jun", 
    "K27.Jul", "K27.Aug", "K27.Sep", "K27.Oct", "K27.Nov", "K27.Dec")

K4K27.normalized <- cbind(K27[,1:7],K4.normalized, K27.normalized, K27[,34:35])

write.csv(K4K27.normalized, 
          file="K4K27_monthmean_ampK4K27over1_maxK27dateorder_normalized.csv", 
          quote=F, row.names=F)

# Heatmap
K27<-read.csv("K4K27_monthmean_ampK4K27over1_maxK27dateorder_normalized.csv",
              header=T,sep=",")
K4 <- t(apply(K27[,8:19],2,rev))
K27 <- t(apply(K27[,20:31],2,rev))

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK4K27over1_maxK27dateorder_normalized.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0,0,0.7)))
image.plot(K27,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K27me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()

pdf("../figs/heatmap/heatmap_K4K27_monthmean_ampK4K27over1_maxK27dateorder_normalized_K4.pdf", 
    width=1.55, height=1.5)
par(mar=c(1,1,1,1))
par(ps=6)
par(cex=1)
cols=colorRamp(c("white",rgb(0.9,0,0)))
image.plot(K4,col=rgb(cols(0:99/99)/255),axes=F,legend.lab="Normalized H3K4me3",
           legend.shrink=0.7,legend.line=1.4,legend.width=0.3,legend.mar=4)
dev.off()



##### Peak shift illustration #####
x<-seq(pi/2,pi*5/2,pi/334)
s<-sin(x)
s6<-sin(x+pi*5/8)
s12<-sin(x+pi*11/8)

dir.create("../figs/others")

pdf("../figs/others/Peak_shift_illustration_maxK4_minK27.pdf",height=0.71,width=1.7)
par(mar=c(1,1.5,0.6,1.6))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)
plot(s12,type="l", col=rgb(0.9,0,0),xlim=c(0,668),ylim=c(-1.2,1.2),
     xlab="",ylab="",ann=F,axes=F)		
par(new=T)
plot(s,type="l", col=rgb(0,0,0.7),xlim=c(0,668),ylim=c(-1.2,1.2),
     xlab="",ylab="",ann=F,axes=F)
day<-c(365,31,28,31,30,31,30,31,31,30,31,30,31)
par(new=T)
plot(NULL, xlim=c(365,730),ylim=c(-0.2,1.2),xlab="",ylab="",ann=F,axes=F)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
       mean(place[9:10]),mean(place[11:12]))
labelx<-c("1","3","5","7","9","11")
axis(side=1,at=place,label=F,tcl=-0.1,pos=-0.2)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=c(-0.083,0.5,1.083),tcl=-0.1,las=2,pos=365,labels=F)
mtext(c("0.0","0.5","1.0"),at=c(-0.083,0.5,1.083),side=2,line=-0.05,las=1)
mtext("Month",side=1,line=-0.05)
mtext("Modification",side=2,line=0.85)
mtext("levels",side=2,line=0.45)
arrows(365,1.2,730,1.2,length=0)
arrows(730,-0.2,730,1.2,length=0)
arrows(365,-0.2,365,1.2,length=0)
arrows(547.5,-0.2,547.5,1.2,length=0,col="black",lwd=0.5)
arrows(479.0625,-0.2,479.0625,1.2,length=0,col="black",lwd=0.5)
arrows(545.5,0.5,481.0625,0.5,length=0.02,lwd=1,code=3)
par(xpd=T)
arrows(513.2812,0.5,513.2812,1.5,length=0,col="black",lwd=0.5)
arrows(513.2812,1.5,528,1.5,length=0,col="black",lwd=0.5)
text(657,1.5,expression(paste(Delta,"(minK27-maxK4)",sep="")))
par(ps=5)
text(584,0.7,"K4",col=rgb(0.9,0,0))
text(659,0.3,"K27",col=rgb(0,0,0.7))
dev.off()


pdf("../figs/others/Peak_shift_illustration_minK4_maxK27.pdf",height=0.71,width=1.7)
par(mar=c(1,1.5,0.6,1.6))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)
plot(s,type="l", col=rgb(0.9,0,0),xlim=c(0,668),ylim=c(-1.2,1.2),
     xlab="",ylab="",ann=F,axes=F)		
par(new=T)
plot(s6,type="l", col=rgb(0,0,0.7),xlim=c(0,668),ylim=c(-1.2,1.2),
     xlab="",ylab="",ann=F,axes=F)
day<-c(365,31,28,31,30,31,30,31,31,30,31,30,31)
par(new=T)
plot(NULL, xlim=c(365,730),ylim=c(-0.2,1.2),xlab="",ylab="",ann=F,axes=F)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
       mean(place[9:10]),mean(place[11:12]))
labelx<-c("1","3","5","7","9","11")
axis(side=1,at=place,label=F,tcl=-0.1,pos=-0.2)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=c(-0.083,0.5,1.083),tcl=-0.1,las=2,pos=365,labels=F)
mtext(c("0.0","0.5","1.0"),at=c(-0.083,0.5,1.083),side=2,line=-0.05,las=1)
mtext("Month",side=1,line=-0.05)
mtext("Modification",side=2,line=0.85)
mtext("levels",side=2,line=0.45)
arrows(365,1.2,730,1.2,length=0)
arrows(730,-0.2,730,1.2,length=0)
arrows(365,-0.2,365,1.2,length=0)
arrows(547.5,-0.2,547.5,1.2,length=0,col="black",lwd=0.5)
arrows(615.9375,-0.2,615.9375,1.2,length=0,col="black",lwd=0.5)
arrows(549.5,0.5,613.9375,0.5,length=0.02,lwd=1,code=3)
par(xpd=T)
arrows(581.7188,0.5,581.7188,1.35,length=0,col="black",lwd=0.5)
text(657,1.5,expression(paste(Delta,"(maxK27-minK4)",sep="")))
par(ps=5)
text(460,0.8,"K4",col=rgb(0.9,0,0))
text(420,0.2,"K27",col=rgb(0,0,0.7))
dev.off()



##### Peak shift between K4 and K27 #####

# Calculation of peak shift
K4.K27 <- read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")
nrow(K4.K27)

shift <- matrix(NA, nrow=nrow(K4.K27), ncol=6)
colnames(shift) <- c("maxK27-minK4","minK27-maxK4","maxK27-minK27",
                     "minK27-maxK27","maxK4-minK4","minK4-maxK4")

for(i in 1:nrow(K4.K27)){
  if(abs(K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) <= 182){
    shift[i,1] <- 
      K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]
  }else if((K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) > 0){
    shift[i,1] <- 
      abs(K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) - 365
  }else if((K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) < 0){
    shift[i,1] <- 
      365 - abs(K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i])
  }else if((K4.K27$K27.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) == 0){
    shift[i,1] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if(abs(K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) <= 182){
    shift[i,2] <- 
      K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]
  }else if((K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) > 0){
    shift[i,2] <- 
      abs(K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) - 365
  }else if((K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) < 0){
    shift[i,2] <- 
      365 - abs(K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i])
  }else if((K4.K27$K27.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) == 0){
    shift[i,2] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$K27.max.date.fromJan1[i] - K4.K27$K27.min.date.fromJan1[i]) > 0){
    shift[i,3] <- 
      K4.K27$K27.max.date.fromJan1[i] - K4.K27$K27.min.date.fromJan1[i]
  }else if((K4.K27$K27.max.date.fromJan1[i] - K4.K27$K27.min.date.fromJan1[i]) < 0){
    shift[i,3] <- 
      365 - abs(K4.K27$K27.max.date.fromJan1[i] - K4.K27$K27.min.date.fromJan1[i])
  }else if((K4.K27$K27.max.date.fromJan1[i] - K4.K27$K27.min.date.fromJan1[i]) == 0){
    shift[i,3] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$K27.min.date.fromJan1[i] - K4.K27$K27.max.date.fromJan1[i]) > 0){
    shift[i,4] <- 
      K4.K27$K27.min.date.fromJan1[i] - K4.K27$K27.max.date.fromJan1[i]
  }else if((K4.K27$K27.min.date.fromJan1[i] - K4.K27$K27.max.date.fromJan1[i]) < 0){
    shift[i,4] <- 
      365 - abs(K4.K27$K27.max.date.fromJan1[i] - K4.K27$K27.min.date.fromJan1[i])
  }else if((K4.K27$K27.min.date.fromJan1[i] - K4.K27$K27.max.date.fromJan1[i]) == 0){
    shift[i,4] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$K4.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) > 0){
    shift[i,5] <- 
      K4.K27$K4.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]
  }else if((K4.K27$K4.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) < 0){
    shift[i,5] <- 
      365 - abs(K4.K27$K4.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i])
  }else if((K4.K27$K4.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i]) == 0){
    shift[i,5] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$K4.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) > 0){
    shift[i,6] <- 
      K4.K27$K4.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]
  }else if((K4.K27$K4.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) < 0){
    shift[i,6] <- 
      365 - abs(K4.K27$K4.max.date.fromJan1[i] - K4.K27$K4.min.date.fromJan1[i])
  }else if((K4.K27$K4.min.date.fromJan1[i] - K4.K27$K4.max.date.fromJan1[i]) == 0){
    shift[i,6] <- 0
  }
}

K4.K27.shift <- cbind(K4.K27,shift)
write.csv(K4.K27.shift, 
          file="K4K27_monthmean_ampK4K27over1_peakshift.csv", quote=F, row.names=F)


# Histogram
K4.K27.shift.over1 <- 
  read.csv("K4K27_monthmean_ampK4K27over1_peakshift.csv",header=T,sep=",")

pdf("../figs/histogram/histogram_maxK27minK4_K4K27_monthmean_ampK4K27over1_peakshift.pdf",
    width=1.4,height=1.3)
par(mar=c(2.2,1.5,2,0.4))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
label <- c("-180","0","180")
hist(K4.K27.shift.over1$maxK27.minK4, las=1, tcl=-0.1, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-200,200), ylim=c(0,120),breaks=12,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0.1,0), at=seq(-180,180,60), labels=F, pos=-10)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,120,60))
mtext(label,side=1,at=c(-200,0,180),line=-0.25)
mtext(expression(paste(Delta,"(maxK27-minK4)",sep="")),side=1,line=0.2)
mtext("No. of genes",side=2,line=0.75)
arrows(0,0,0,120,col="orange",length=0)
arrows(0,112,50,112,length=0.03,col="orange")
par(ps=5)
text(110,112,"Delayed",col="orange")
text(112,90,"K27",col="orange")
arrows(0,112,-50,112,length=0.03,col="orange")
text(-122,112,"Advanced",col="orange")
text(-112,90,"K27",col="orange")
dev.off()


pdf("../figs/histogram/histogram_minK27maxK4_K4K27_monthmean_ampK4K27over1_peakshift.pdf",
    width=1.4,height=1.3)
par(mar=c(2.2,1.5,2,0.4))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(K4.K27.shift.over1$minK27.maxK4, las=1, tcl=-0.1, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-200,200), ylim=c(0,120),breaks=24,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0.1,0), at=seq(-180,180,60), labels=F, pos=-10)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,120,60))
mtext(label,side=1,at=c(-200,0,180),line=-0.25)
mtext(expression(paste(Delta,"(minK27-maxK4)",sep="")),side=1,line=0.2)
mtext("No. of genes",side=2,line=0.75)
arrows(0,0,0,120,col="orange",length=0)
arrows(0,112,50,112,length=0.03,col="orange")
par(ps=5)
text(110,112,"Delayed",col="orange")
text(112,90,"K27",col="orange")
arrows(0,112,-50,112,length=0.03,col="orange")
text(-122,112,"Advanced",col="orange")
text(-112,90,"K27",col="orange")
dev.off()



##### Peak shift between K4 and K27 for cosinor model #####

# Calculation of peak shift
K4.K27 <- read.csv("K4K27_monthmean_ampK4K27over1.csv",header=T,sep=",")
cosinor <- read.csv("K4K27_cosinor.csv",header=T,sep=",")
K4.K27 <- merge(K4.K27,cosinor, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type",
             "Symbol","Alias","Note"),all=F,sort=F)
nrow(K4.K27)
nrow(cosinor)
head(K4.K27)

shift <- matrix(NA, nrow=nrow(K4.K27), ncol=6)
colnames(shift) <- c("maxK27-minK4","minK27-maxK4","maxK27-minK27",
                     "minK27-maxK27","maxK4-minK4","minK4-maxK4")

for(i in 1:nrow(K4.K27)){
  if(abs(K4.K27$peak.K27[i] - K4.K27$trough.K4[i]) <= 182){
    shift[i,1] <- 
      K4.K27$peak.K27[i] - K4.K27$trough.K4[i]
  }else if((K4.K27$peak.K27[i] - K4.K27$trough.K4[i]) > 0){
    shift[i,1] <- 
      abs(K4.K27$peak.K27[i] - K4.K27$trough.K4[i]) - 365
  }else if((K4.K27$peak.K27[i] - K4.K27$trough.K4[i]) < 0){
    shift[i,1] <- 
      365 - abs(K4.K27$peak.K27[i] - K4.K27$trough.K4[i])
  }else if((K4.K27$peak.K27[i] - K4.K27$trough.K4[i]) == 0){
    shift[i,1] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if(abs(K4.K27$trough.K27[i] - K4.K27$peak.K4[i]) <= 182){
    shift[i,2] <- 
      K4.K27$trough.K27[i] - K4.K27$peak.K4[i]
  }else if((K4.K27$trough.K27[i] - K4.K27$peak.K4[i]) > 0){
    shift[i,2] <- 
      abs(K4.K27$trough.K27[i] - K4.K27$peak.K4[i]) - 365
  }else if((K4.K27$trough.K27[i] - K4.K27$peak.K4[i]) < 0){
    shift[i,2] <- 
      365 - abs(K4.K27$trough.K27[i] - K4.K27$peak.K4[i])
  }else if((K4.K27$trough.K27[i] - K4.K27$peak.K4[i]) == 0){
    shift[i,2] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$peak.K27[i] - K4.K27$trough.K27[i]) > 0){
    shift[i,3] <- 
      K4.K27$peak.K27[i] - K4.K27$trough.K27[i]
  }else if((K4.K27$peak.K27[i] - K4.K27$trough.K27[i]) < 0){
    shift[i,3] <- 
      365 - abs(K4.K27$peak.K27[i] - K4.K27$trough.K27[i])
  }else if((K4.K27$peak.K27[i] - K4.K27$trough.K27[i]) == 0){
    shift[i,3] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$trough.K27[i] - K4.K27$peak.K27[i]) > 0){
    shift[i,4] <- 
      K4.K27$trough.K27[i] - K4.K27$peak.K27[i]
  }else if((K4.K27$trough.K27[i] - K4.K27$peak.K27[i]) < 0){
    shift[i,4] <- 
      365 - abs(K4.K27$peak.K27[i] - K4.K27$trough.K27[i])
  }else if((K4.K27$trough.K27[i] - K4.K27$peak.K27[i]) == 0){
    shift[i,4] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$peak.K4[i] - K4.K27$trough.K4[i]) > 0){
    shift[i,5] <- 
      K4.K27$peak.K4[i] - K4.K27$trough.K4[i]
  }else if((K4.K27$peak.K4[i] - K4.K27$trough.K4[i]) < 0){
    shift[i,5] <- 
      365 - abs(K4.K27$peak.K4[i] - K4.K27$trough.K4[i])
  }else if((K4.K27$peak.K4[i] - K4.K27$trough.K4[i]) == 0){
    shift[i,5] <- 0
  }
}

for(i in 1:nrow(K4.K27)){
  if((K4.K27$trough.K4[i] - K4.K27$peak.K4[i]) > 0){
    shift[i,6] <- 
      K4.K27$trough.K4[i] - K4.K27$peak.K4[i]
  }else if((K4.K27$trough.K4[i] - K4.K27$peak.K4[i]) < 0){
    shift[i,6] <- 
      365 - abs(K4.K27$peak.K4[i] - K4.K27$trough.K4[i])
  }else if((K4.K27$trough.K4[i] - K4.K27$peak.K4[i]) == 0){
    shift[i,6] <- 0
  }
}

K4.K27.shift <- cbind(K4.K27,shift)
write.csv(K4.K27.shift, 
          file="K4K27_monthmean_ampK4K27over1_cosinor_peakshift.csv", quote=F, row.names=F)


# Histogram
K4.K27.shift.over1 <- 
  read.csv("K4K27_monthmean_ampK4K27over1_cosinor_peakshift.csv",header=T,sep=",")

pdf("../figs/histogram/histogram_maxK27minK4_K4K27_monthmean_ampK4K27over1_cosinor_peakshift.pdf",
    width=1.4,height=1.3)
par(mar=c(2.2,1.5,2,0.4))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
label <- c("-180","0","180")
hist(K4.K27.shift.over1$maxK27.minK4, las=1, tcl=-0.1, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-200,200), ylim=c(0,160),breaks=12,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0.1,0), at=seq(-180,180,60), labels=F, pos=-10*4/3)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,160,80))
mtext(label,side=1,at=c(-200,0,180),line=-0.25)
mtext(expression(paste(Delta,"(peak K27 - trough K4)",sep="")),side=1,line=0.2)
mtext("No. of genes",side=2,line=0.75)
arrows(0,0,0,160,col="orange",length=0)
arrows(0,112*4/3,50,112*4/3,length=0.03,col="orange")
par(ps=5)
text(110,112*4/3,"Delayed",col="orange")
text(112,90*4/3,"K27",col="orange")
arrows(0,112*4/3,-50,112*4/3,length=0.03,col="orange")
text(-122,112*4/3,"Advanced",col="orange")
text(-112,90*4/3,"K27",col="orange")
dev.off()


pdf("../figs/histogram/histogram_minK27maxK4_K4K27_monthmean_ampK4K27over1_cosinor_peakshift.pdf",
    width=1.4,height=1.3)
par(mar=c(2.2,1.5,2,0.4))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(K4.K27.shift.over1$minK27.maxK4, las=1, tcl=-0.1, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-200,200), ylim=c(0,160),breaks=12,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0.1,0), at=seq(-180,180,60), labels=F, pos=-10*4/3)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,160,80))
mtext(label,side=1,at=c(-200,0,180),line=-0.25)
mtext(expression(paste(Delta,"(trough K27 - peak K4)",sep="")),side=1,line=0.2)
mtext("No. of genes",side=2,line=0.75)
arrows(0,0,0,160,col="orange",length=0)
arrows(0,112*4/3,50,112*4/3,length=0.03,col="orange")
par(ps=5)
text(110,112*4/3,"Delayed",col="orange")
text(112,90*4/3,"K27",col="orange")
arrows(0,112*4/3,-50,112*4/3,length=0.03,col="orange")
text(-122,112*4/3,"Advanced",col="orange")
text(-112,90*4/3,"K27",col="orange")
dev.off()



##### Preparation for spline to calculate the area of Lissajous curves #####
K4.K27.duplicate <- 
  read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv",header=T,sep=",")
K4 <- 
  cbind(K4.K27.duplicate[,11:22],K4.K27.duplicate[,11:22],K4.K27.duplicate[,11:22],
        K4.K27.duplicate[,41:52],K4.K27.duplicate[,41:52],K4.K27.duplicate[,41:52])
K27 <- 
  cbind(K4.K27.duplicate[,26:37],K4.K27.duplicate[,26:37],K4.K27.duplicate[,26:37],
        K4.K27.duplicate[,56:67],K4.K27.duplicate[,56:67],K4.K27.duplicate[,56:67])
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
nrow(K4)
nrow(K27)

date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)



##### Estimation of max and min of plot region #####
maxK4 <- rep(NA,length=nrow(K4))
minK4 <- rep(NA,length=nrow(K4))
maxK27 <- rep(NA,length=nrow(K4))
minK27 <- rep(NA,length=nrow(K4))

for(i in 1:nrow(K4)){
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  maxK4[i] <- max(predK4$y)
  minK4[i] <- min(predK4$y)
  maxK27[i] <- max(predK27$y)
  minK27[i] <- min(predK27$y)
}

max(maxK4)
min(minK4)
max(maxK27)
min(minK27)



##### Drawing of Lissajous curves to calculate the area #####
dir.create("../figs/Lissajous_naked")
for(i in 1:nrow(K4)){
  
  pdf(paste("../figs/Lissajous_naked/Lissajous_maxK4K27over2_ampK4K27over1_", 
            i, ".pdf", sep=""),width=2,height=2.6)
  
  par(xpd=NA)
  
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  plot(predK4$y,predK27$y,xlab="",ylab="",main="",axes=F,ann=F,
       type="l",lwd=0.5,xlim=c(0,9),ylim=c(0,9))
  
  dev.off()
}



##### Calculation of the area of Lissajous curves by Adobe illustrater #####

# Preparation of the csv file with the format below
# Max.Area.millimeters, row
# 7.666, 1
# 1.129, 10
# 2.969, 100
# ...........
# Save as "Lissajous_maxK4K27over2_ampK4K27over1_area.csv" in "data" folder



##### Addition of the area to the last column of the data file and ordering #####
max.area <- 
  read.csv("Lissajous_maxK4K27over2_ampK4K27over1_area.csv",sep=",",header=T)
sortlist <- order(max.area$row)
max.area <- max.area[sortlist,]

K4.K27 <- 
  read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv",header=T,sep=",")
K4.K27.max.area <- cbind(K4.K27, max.area)

sortlist <- order(K4.K27.max.area$Max.Area.millimeters, decreasing=T)
K4.K27.max.area.order <- K4.K27.max.area[sortlist,]
K4.K27.max.area.order <- K4.K27.max.area.order[,-ncol(K4.K27.max.area.order)]

write.csv(K4.K27.max.area.order, 
          file="K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.csv",
          quote=F, row.names=F)

K4.K27.max.area.order2 <- K4.K27.max.area.order[,c(1:7,23:25,68:78)]
write.csv(K4.K27.max.area.order2, 
          file="K4K27_stat_maxK4K27over2_ampK4K27over1_area.csv",
          quote=F, row.names=F)


# Overwriting the file
max.area <- read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.csv",sep=",",header=T)
K4.K27 <- read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv",sep=",",header=T)
max.area.K4.K27 <- merge(max.area, K4.K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                "Locus_type","Symbol","Alias","Note"),
                         all=F,sort=F)
max.area.K4.K27 <- max.area.K4.K27[,1:ncol(max.area)]
colnames(max.area.K4.K27) <- gsub(".x", "", colnames(max.area.K4.K27))
colnames(max.area.K4.K27) <- gsub("mK", "maxK", colnames(max.area.K4.K27))
colnames(max.area.K4.K27) <- gsub("M.A", "Max.A", colnames(max.area.K4.K27))
nrow(max.area)
nrow(K4.K27)
nrow(max.area.K4.K27)
write.csv(max.area.K4.K27, 
          file="K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.csv",
          quote=F, row.names=F)

max.area.K4.K27.2 <- max.area.K4.K27[,c(1:7,23:25,68:78)]
write.csv(max.area.K4.K27.2, 
          file="K4K27_stat_maxK4K27over2_ampK4K27over1_area.csv",
          quote=F, row.names=F)



##### Lissajous curves of K4-K27 for FLC #####
K4.K27 <- 
  read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.csv",
           header=T,sep=",")

K4 <- cbind(K4.K27[,11:22],K4.K27[,11:22],K4.K27[,11:22],
            K4.K27[,41:52],K4.K27[,41:52],K4.K27[,41:52])
K27 <- cbind(K4.K27[,26:37],K4.K27[,26:37],K4.K27[,26:37],
             K4.K27[,56:67],K4.K27[,56:67],K4.K27[,56:67])
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)

date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27[,1:7]

dir.create("../figs/Lissajous_curve")

pdf("../figs/Lissajous_curve/Lissajous_FLC_K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.pdf",
    width=1.1,height=1.1)
for(i in 1){
  par(ps=6)
  par(mar=c(2,2,1,1))
  par(mgp=c(0,0.3,0))
  plot(K4[,c(1:12,37:48)], K27[,c(1:12,37:48)], pch=".", col="gray",
       las=1, tcl=-0.2,xlim=c(0,9),ylim=c(0,9), main="",axes=F,ann=F)
  par(new=T)
  plot(K4[i,c(1:12,37:48)], K27[i,c(1:12,37:48)],  pch=".",las=1,
       tcl=-0.2,xlim=c(0,9),ylim=c(0,9),axes=F, ann=F,xlab="",ylab="",type="p")
  
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365+6,730+6,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  par(new=T)
  plot(predK4$y,predK27$y,xlab="",ylab="",main="",axes=F,ann=F,
       type="l",lwd=0.5,xlim=c(0,9),ylim=c(0,9), col="red")
  label <- c("0","3","6","9")
  box()
  axis(1, at=c(0,3,6,9), tcl=-0.1,labels=F)
  axis(2, at=c(0,3,6,9), tcl=-0.1,labels=F)
  mtext(label,side=1,at=c(0,3,6,9),line=-0.33)
  mtext(label,side=2,at=c(0,3,6,9),line=0.15,las=1)
  mtext("H3K27me3",side=2,line=0.8)
  mtext(expression(paste("[",log[2],"(rpkm)]")),side=2,line=0.3)
  mtext("H3K4me3",side=1,line=0.1)
  mtext(expression(paste("[",log[2],"(rpkm)]")),side=1,line=0.6)
  if(is.na(title[i,5])){mtext(paste(title[i,1], " : ", "N.A.", sep=""),
              side=3,line=0.1,font=3,at=3)
  }else{mtext(paste(title[i,1], " : ", "Ahg", title[i,5], sep=""),
                               side=3,line=0.1,font=3,at=3)
  }
}
dev.off()



##### Lissajous curves of K4-K27 for other genes #####
for(i in 2:4){
pdf(paste("../figs/Lissajous_curve/Lissajous_",title[i,1],"_","Ahg",title[i,5], "_K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.pdf",sep=""),
    width=1.3,height=1.3)
  par(ps=6)
  par(mar=c(2,2,1,1))
  par(mgp=c(0,0.3,0))
  plot(K4[,c(1:12,37:48)], K27[,c(1:12,37:48)], pch=".", col="gray",
       las=1, tcl=-0.2,xlim=c(0,9),ylim=c(0,9), main="",axes=F,ann=F)
  par(new=T)
  plot(K4[i,c(1:12,37:48)], K27[i,c(1:12,37:48)],  pch=".",las=1, 
       tcl=-0.2,xlim=c(0,9),ylim=c(0,9),axes=F, ann=F,xlab="",ylab="",type="p")
  
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365+6,730+6,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)
  par(new=T)
  plot(predK4$y,predK27$y,xlab="",ylab="",main="",axes=F,ann=F,
       type="l",lwd=0.5,xlim=c(0,9),ylim=c(0,9), col="red")
  label <- c("0","3","6","9")
  box()
  axis(1, at=c(0,3,6,9), tcl=-0.1, mgp=c(0,-0.1,0),labels=F)
  mtext(label,side=1,at=c(0,3,6,9),line=-0.3)
  axis(2, at=c(0,3,6,9), tcl=-0.1, mgp=c(0,-0.1,0),labels=F)
  mtext(label,side=2,at=c(0,3,6,9),line=0.15,las=1)
  mtext("H3K27me3",side=2,line=0.9)
  mtext(expression(paste("[",log[2],"(rpkm)]")),side=2,line=0.4)
  mtext("H3K4me3",side=1,line=0.2)
  mtext(expression(paste("[",log[2],"(rpkm)]")),side=1,line=0.7)
  if(is.na(title[i,5])){mtext(paste(title[i,1], " : ", "N.A.", sep=""),
              side=3,line=0.1,font=3,at=4.5)
  }else{mtext(paste(title[i,1], " : ", "Ahg", title[i,5], sep=""),
                               side=3,line=0.1,font=3,at=4.5)
  }
dev.off()
}



##### Determination of the direction of rotation #####
predK4.date <- matrix(NA,nrow=nrow(K4),ncol=12)
predK27.date <- matrix(NA,nrow=nrow(K4),ncol=12)
sampling.date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)

dir.create("../figs/Lissajous_curve/direction")

# All dual seasonality genes
pdf("../figs/Lissajous_curve/direction/Lissajous_Direction_K4K27_data_stat_maxK4K27over2_ampK4K27over1_area.pdf",
    width=1.2,height=1.15)

for(i in 1:nrow(K4)){
  
  for(j in 1:12){
    spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
    spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
    length <- 1000
    x <- seq(365,730,length=length)
    predK4 <- predict(spK4,x)
    predK27 <- predict(spK27,x)
    predK4$y <- (predK4$y-min(predK4$y))/(max(predK4$y)-min(predK4$y))
    predK27$y <- (predK27$y-min(predK27$y))/(max(predK27$y)-min(predK27$y))
    
    predK4.date[i,j] <- 
      predK4$y[which.min(abs(predK4$x - (365+sampling.date[j])))]
    predK27.date[i,j] <- 
      predK27$y[which.min(abs(predK27$x - (365+sampling.date[j])))]
  }
  
  par(ps=6)
  par(mar=c(1.5,1.7,1,1))
  par(mgp=c(0,0.3,0))
  
  plot(predK4.date[i,1], predK27.date[i,1], pch=1, las=1, 
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="",cex=0.3,axes=F,ann=F, lwd=0.3)
  par(new=T)
  plot(predK4.date[i,2:4], predK27.date[i,2:4], pch=16, las=1,
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.3, lwd=0.3)
  par(new=T)
  plot(predK4.date[i,5:8], predK27.date[i,5:8], pch=0, las=1,
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.3, lwd=0.3)
  par(new=T)
  plot(predK4.date[i,9:11], predK27.date[i,9:11], pch=15, las=1, 
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.3, lwd=0.3)
  par(new=T)
  plot(predK4.date[i,12], predK27.date[i,12], pch=1, las=1, 
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.3, lwd=0.3)
  par(new=T)
  plot(predK4.date[i,1:12], predK27.date[i,1:12], type="l", 
       lwd=0.3, las=1, tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.3)
  par(new=T)
  plot(predK4.date[i,c(12,1)], predK27.date[i,c(12,1)],type="l",
       lty="dotted", lwd=0.3,las=1, tcl=-0.2, xlab="",ylab="",xlim=c(0,1),
       ylim=c(0,1), main="", axes=F,ann=F,cex=0.3)
  
  par(xpd=T)
  textxy(predK4.date[i,1:12], predK27.date[i,1:12], 1:12, cex=0.8)
 
  box()
  axis(1, at=c(0,0.5,1), tcl=-0.1, labels=F)
  axis(2, at=c(0,0.5,1), tcl=-0.1, labels=F)
  mtext(c("0.0","0.5","1.0"),at=c(0,0.5,1),side=1,line=-0.35)
  mtext(c("0.0","0.5","1.0"),at=c(0,0.5,1),side=2,line=0.15,las=1)
  mtext("H3K27me3",side=2,line=0.7)
  mtext("H3K4me3",side=1,line=0.1)
  
  if(is.na(title[i,5])){mtext(paste(title[i,1], " : ", "N.A.", sep=""),
              side=3,line=0.1,font=3)
  }else{mtext(paste(title[i,1], " : ", "Ahg", title[i,5], sep=""),
                               side=3,line=0.1,font=3)
  }
  
}
dev.off()



##### Area top 20 genes #####
pdf("../figs/Lissajous_curve/direction/Lissajous_Direction_K4K27_data_stat_maxK4K27over2_ampK4K27over1_area_top20genes.pdf",
    width=5.88,height=7.7)
# width=4.2,height=5.5,cex=0.3,lwd=0.3 for small figures
  set.panel(5,4,relax=T)
  par(ps=6)
  par(oma=c(0,0,0,0))
  par(mar=c(1.2,1.5,1,0.4))
  par(mgp=c(0,0.3,0))
  par(xpd=T)
  par(cex=T)

for(i in 1:20){
  
  for(j in 1:12){
    spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
    spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
    length <- 1000
    x <- seq(365,730,length=length)
    predK4 <- predict(spK4,x)
    predK27 <- predict(spK27,x)
    predK4$y <- (predK4$y-min(predK4$y))/(max(predK4$y)-min(predK4$y))
    predK27$y <- (predK27$y-min(predK27$y))/(max(predK27$y)-min(predK27$y))
    
    predK4.date[i,j] <- 
      predK4$y[which.min(abs(predK4$x - (365+sampling.date[j])))]
    predK27.date[i,j] <- 
      predK27$y[which.min(abs(predK27$x - (365+sampling.date[j])))]
  }
  
  plot(predK4.date[i,1], predK27.date[i,1], pch=1, las=1, 
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="",cex=0.8,axes=F,ann=F, lwd=0.6)
  par(new=T)
  plot(predK4.date[i,2:4], predK27.date[i,2:4], pch=16, las=1,
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.8, lwd=0.6)
  par(new=T)
  plot(predK4.date[i,5:8], predK27.date[i,5:8], pch=0, las=1,
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.8, lwd=0.6)
  par(new=T)
  plot(predK4.date[i,9:11], predK27.date[i,9:11], pch=15, las=1, 
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.8, lwd=0.6)
  par(new=T)
  plot(predK4.date[i,12], predK27.date[i,12], pch=1, las=1, 
       tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.8, lwd=0.6)
  par(new=T)
  plot(predK4.date[i,1:12], predK27.date[i,1:12], type="l", 
       lwd=0.6, las=1, tcl=-0.2, xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), 
       main="", axes=F,ann=F,cex=0.8)
  par(new=T)
  plot(predK4.date[i,c(12,1)], predK27.date[i,c(12,1)],type="l",
       lty="dotted", lwd=0.6,las=1, tcl=-0.2, xlab="",ylab="",xlim=c(0,1),
       ylim=c(0,1), main="", axes=F,ann=F,cex=0.8)
  
  textxy(predK4.date[i,1:12], predK27.date[i,1:12], 1:12, cex=1, offset=1.2)	
  
  box()
  axis(1, at=c(0,0.5,1), tcl=-0.1, labels=F)
  axis(2, at=c(0,0.5,1), tcl=-0.1, labels=F)
  mtext(c("0.0","0.5","1.0"),at=c(0,0.5,1),side=1,line=-0.35)
  mtext(c("0.0","0.5","1.0"),at=c(0,0.5,1),side=2,line=0.15,las=1)
  mtext("H3K27me3",side=2,line=0.7)
  mtext("H3K4me3",side=1,line=0.1)
  
  if(is.na(title[i,5])){mtext(paste(title[i,1], " : ", "N.A.", sep=""),
              side=3,line=0.1,font=3)
  }else{mtext(paste(title[i,1], " : ", "Ahg", title[i,5], sep=""),
                               side=3,line=0.1,font=3)
  }
  
}
set.panel()
dev.off()



##### Calculation of relative area #####
K4.K27.cor <- 
  read.csv("K4K27_ampK4K27over1_cor.csv",header=T,sep=",")
K4.K27.area <- 
  read.csv("K4K27_stat_maxK4K27over2_ampK4K27over1_area.csv",header=T,sep=",")
K4.K27.area.cor <- 
  merge(K4.K27.area, K4.K27.cor, 
        by=c("Ahal_ID","Araport11_ID",
             "Reciprocal_flag","Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)

Relative.area <- 
  K4.K27.area.cor$Max.Area.millimeters / max(K4.K27.area.cor$Max.Area.millimeters)
K4.K27.area.cor.relative <- 
  cbind(K4.K27.area.cor[,-(ncol(K4.K27.area.cor)-1)], Relative.area)

colnames(K4.K27.area.cor.relative)[8:10] <- c("chr","start","end")

write.csv(K4.K27.area.cor.relative, 
          file="K4K27_stat_maxK4K27over2_ampK4K27over1_cor_relativearea.csv",
          quote=F, row.names=F)



###### Histogram of relative area #####
# In advance, add "Direction" to the last column of the file, ..._relativearea.csv,
# and save it as ..._relativearea_dir.csv

area <- 
  read.csv("K4K27_stat_maxK4K27over2_ampK4K27over1_cor_relativearea_dir.csv",
           header=T,sep=",")

area.mod <- rep(NA, length=nrow(area))
for(i in 1:nrow(area)){
  if(area$Relative.area[i] <= 0.5){
    area.mod[i] <- area$Relative.area[i]
  }else{area.mod[i] <- 0.5001}
}

area.over0.15 <- area.mod[1:60]

area.clock <- rep(NA, length=60)
for(i in 1:60){
  if(area$Direction[i]=="Clock"){
    if(area$Relative.area[i] <= 0.5){
      area.clock[i] <- area$Relative.area[i]
    }else{area.clock[i] <- 0.5001}
  }
}

pdf("../figs/histogram/histogram_K4K27_stat_maxK4K27over2_ampK4K27over1_cor_relativearea_dir.pdf",
    width=1, height=0.8)
par(mar=c(1,1.5,0.5,0.5))
par(ps=6)
par(xpd=T)
xlab <- c("0.0","0.2","0.4")
ylab <- c("0","60","120")
hist(area.mod, main="", xlim=c(0, 0.55), ylim=c(0,120), las=1, axes=F,
     xlab="",ylab="",col="white", breaks=12)
hist(area.over0.15, main="", xlim=c(0, 0.55), ylim=c(0,120), las=1, axes=F,
     xlab="",ylab="",col="gray", breaks=6, add=T)
hist(area.clock, main="", xlim=c(0, 0.55), ylim=c(0,120), las=1, axes=F,
     xlab="",ylab="",col="darkorange", breaks=6, add=T)
axis(side=1,at=seq(0,0.5,0.1),tcl=-0.1,las=1,mgp=c(0,0,0),labels=F)
axis(side=2,at=seq(0,120,60),tcl=-0.1,las=1,mgp=c(0,0.3,0),labels=F)
rect(0,0,0.55,30,border="red")
mtext(xlab,side=1,at=c(0,0.2,0.4),line=-0.35,las=1)
mtext(ylab,side=2,at=seq(0,120,60),line=0.15,las=1)
dev.off()


pdf("../figs/histogram/histogram_magnif_K4K27_stat_maxK4K27over2_ampK4K27over1_cor_relativearea_dir.pdf",
    width=1.5, height=1)
par(mar=c(1.5,1.3,0.5,2.2))
par(ps=6)
hist(area.mod, main="", xlim=c(0, 0.55), ylim=c(0,30), las=1, axes=F,
     xlab="",ylab="",col="white", breaks=12)
hist(area.over0.15, main="", xlim=c(0, 0.55), ylim=c(0,30), las=1, axes=F,
     xlab="",ylab="",col="gray", breaks=6, add=T)
hist(area.clock, main="", xlim=c(0, 0.55), ylim=c(0,30), las=1, axes=F,
     xlab="",ylab="",col="darkorange", breaks=6, add=T)
xlab <- c("0.0","0.2","0.4")
ylab <- c("0","10","20","30")
axis(side=1,at=seq(0,0.5,0.1),tcl=-0.1,las=1,mgp=c(0,0,0),labels=F)
axis(side=2,at=seq(0,30,10),tcl=-0.1,las=1,mgp=c(0,0.3,0),labels=F)
mtext(xlab,side=1,at=c(0,0.2,0.4),line=-0.35,las=1)
mtext(ylab,side=2,at=seq(0,30,10),line=0.15,las=1)
mtext("Relative area", side=1, line=0.1)
mtext("No. of genes", side=2, line=0.6)
par(xpd=T)
arrows(0.525,4,0.525,9,length=0,col="black",lwd=0.5)
arrows(0.525,9,0.585,9,length=0,col="black",lwd=0.5)
text(0.7,8.5,expression(italic("AhgFLC")))
text(0.71,3.3,expression(italic("AhgVIN3")))
text(0.724,-1.9,expression(italic("AhgMAF1")))
text(0.739,-7.1,expression(italic("AhgRAB18")))

dev.off()



##### Merge files #####
K4K27.area<-
  read.csv("K4K27_stat_maxK4K27over2_ampK4K27over1_cor_relativearea_dir.csv",
           header=T,sep=",")
K4K27mean <- read.csv("K4K27_monthmean.csv",header=T,sep=",")

K4K27.area.mean <- 
  merge(K4K27.area,K4K27mean, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
             "Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)
K4K27.area.mean <- K4K27.area.mean[,c(1:10,24:ncol(K4K27.area.mean),11:23)]

write.csv(K4K27.area.mean, 
          file="K4K27_monthmean_ampK4K27over1_stat_area.csv", 
          quote=F, row.names=F)
