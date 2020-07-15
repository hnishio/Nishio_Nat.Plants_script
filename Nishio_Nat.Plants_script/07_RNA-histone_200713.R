
library(exactRankTests)
library(lawstat)

source("functions/GOanalysis_functions.R")



##### Preparation for spline #####
setwd("data")

load("20181120_bowtie2_rpkm.Rdata")

genes.log2rpkm <- log2(rpkm+1)
genes.log2rpkm2 <- cbind(genes.log2rpkm,genes.log2rpkm,genes.log2rpkm)
nrow(genes.log2rpkm2)

attribute <- read.table("140519_SampleAttribute.txt", header=T,sep="\t")
attribute.date <- attribute[1:490,2:5]
times <- as.POSIXct(paste(attribute.date[,'year'], 
                          attribute.date[,'month'], 
                          attribute.date[,'day'],
                          attribute.date[,'hour'],
                          sep='/'), 
                    format='%Y/%m/%d/%H')

days <- as.numeric(times - as.POSIXct("2011-06-30 12:00:00"))
days2 <- c(days,730+days,1460+days)


# Calculation of max, min, mean, sd, amp
maxRNA <- rep(NA,length=nrow(genes.log2rpkm2))
minRNA <- rep(NA,length=nrow(genes.log2rpkm2))
meanRNA <- rep(NA,length=nrow(genes.log2rpkm2))
sdRNA <- rep(NA,length=nrow(genes.log2rpkm2))
ampRNA <- rep(NA,length=nrow(genes.log2rpkm2))

for(i in 1:nrow(genes.log2rpkm2)){
spRNA <- smooth.spline(days2,genes.log2rpkm2[i,],spar=0.3)
length <- 10000
x <- seq(730,1460,length=length)
predRNA <- predict(spRNA,x)

#plot(730+days, genes.log2rpkm[19190,])
#lines(predRNA, col=2)

maxRNA[i] <- max(predRNA$y)
minRNA[i] <- min(predRNA$y)
meanRNA[i] <- mean(predRNA$y)
sdRNA[i] <- sd(predRNA$y)
ampRNA[i] <- max(predRNA$y)-min(predRNA$y)
}

Ahal_ID <- rep(NA,length=nrow(genes.log2rpkm2))
list <- strsplit(row.names(genes.log2rpkm2), ".t")
for(i in 1:length(Ahal_ID)){
  Ahal_ID[i] <- list[[i]][1]
}

RNA.stat <- as.data.frame(cbind(Ahal_ID,maxRNA,minRNA,meanRNA,sdRNA,ampRNA))
write.csv(RNA.stat, file="RNA_stat.csv", quote=F, row.names=F)



##### Calculation of monthly mean #####
RNA.month <- matrix(NA,nrow=nrow(genes.log2rpkm2),ncol=12)
colnames(RNA.month) <- 
  c("RNA.Jul", "RNA.Aug", "RNA.Sep", "RNA.Oct", "RNA.Nov", "RNA.Dec", 
    "RNA.Jan", "RNA.Feb", "RNA.Mar", "RNA.Apr", "RNA.May", "RNA.Jun")
RNA.month <- as.data.frame(RNA.month)

RNA.max.date.fromJan1 <- rep(NA,length=nrow(genes.log2rpkm2))
RNA.min.date.fromJan1 <- rep(NA,length=nrow(genes.log2rpkm2))

for(i in 1:nrow(RNA.month)){
  spRNA <- smooth.spline(days2,genes.log2rpkm2[i,],spar=0.3)
  length <- 10000
  x <- seq(730,1460,length=length)
  predRNA <- predict(spRNA,x)
  
  # monthly mean
  RNA.month$RNA.Jul[i] <- mean(predRNA$y[(730<predRNA$x & predRNA$x<=730+31) | 
                                         (730+365<predRNA$x & predRNA$x<=730+365+31)])
  RNA.month$RNA.Aug[i] <- mean(predRNA$y[(730+31<predRNA$x & predRNA$x<=730+62) |
                                         (730+365+31<predRNA$x & predRNA$x<=730+365+62)])
  RNA.month$RNA.Sep[i] <- mean(predRNA$y[(730+62<predRNA$x & predRNA$x<=730+92) |
                                         (730+365+62<predRNA$x & predRNA$x<=730+365+92)])
  RNA.month$RNA.Oct[i] <- mean(predRNA$y[(730+92<predRNA$x & predRNA$x<=730+123) |
                                         (730+365+92<predRNA$x & predRNA$x<=730+365+123)])
  RNA.month$RNA.Nov[i] <- mean(predRNA$y[(730+123<predRNA$x & predRNA$x<=730+153) |
                                         (730+365+123<predRNA$x & predRNA$x<=730+365+153)])
  RNA.month$RNA.Dec[i] <- mean(predRNA$y[(730+153<predRNA$x & predRNA$x<=730+184) |
                                         (730+365+153<predRNA$x & predRNA$x<=730+365+184)])
  RNA.month$RNA.Jan[i] <- mean(predRNA$y[(730+184<predRNA$x & predRNA$x<=730+215) |
                                         (730+365+184<predRNA$x & predRNA$x<=730+365+215)])
  RNA.month$RNA.Feb[i] <- mean(predRNA$y[(730+215<predRNA$x & predRNA$x<=730+243) |
                                         (730+365+215<predRNA$x & predRNA$x<=730+365+243)])
  RNA.month$RNA.Mar[i] <- mean(predRNA$y[(730+243<predRNA$x & predRNA$x<=730+274) |
                                         (730+365+243<predRNA$x & predRNA$x<=730+365+274)])
  RNA.month$RNA.Apr[i] <- mean(predRNA$y[(730+274<predRNA$x & predRNA$x<=730+304) |
                                         (730+365+274<predRNA$x & predRNA$x<=730+365+304)])
  RNA.month$RNA.May[i] <- mean(predRNA$y[(730+304<predRNA$x & predRNA$x<=730+335) |
                                         (730+365+304<predRNA$x & predRNA$x<=730+365+335)])
  RNA.month$RNA.Jun[i] <- mean(predRNA$y[(730+335<predRNA$x & predRNA$x<=730+365) |
                                         (730+365+335<predRNA$x & predRNA$x<=730+365+365)])
  
  # Days from 1 Jan.
  spline.2y <- cbind(predRNA$y[730<=predRNA$x & predRNA$x<=730+365],
                     predRNA$y[730+365<=predRNA$x & predRNA$x<=730+730])
  spline.2y.mean <- apply(spline.2y, 1, mean)
  mean.list <- list(predRNA$x[730<=predRNA$x & predRNA$x<=730+365],
                    spline.2y.mean)
  
  max.day <- mean.list[[1]][mean.list[[2]] == max(mean.list[[2]])]
  min.day <- mean.list[[1]][mean.list[[2]] == min(mean.list[[2]])]
  
  if(max(mean.list[[2]])>0){
    if((730+184 < max.day) && (max.day <= 730+365)){
      RNA.max.date.fromJan1[i] <- max.day - 730 - 184
    }else{RNA.max.date.fromJan1[i] <- max.day + 365 - 730 - 184}
    
    if((730+184 < min.day) && (min.day <= 730+365)){
      RNA.min.date.fromJan1[i] <- min.day - 730 - 184
    }else{RNA.min.date.fromJan1[i] <- min.day + 365 - 730 - 184}
  }else{}
}

RNA.mat <- cbind(Ahal_ID, RNA.month, 
                 RNA.max.date.fromJan1, RNA.min.date.fromJan1)
RNA.mat2 <- RNA.mat[,c(1,8:13,2:7,14:ncol(RNA.mat))]
 
write.csv(RNA.mat2, file="RNA_monthmean.csv", quote=F, row.names=F)



##### Extraction of genes with max > 2 #####
RNA <- read.csv("RNA_stat.csv",header=T,sep=",")
RNA.over2 <- subset(RNA,maxRNA>2)
nrow(RNA.over2)
write.csv(RNA.over2, 
          file="RNA_stat_maxRNAover2.csv", 
          quote=F, row.names=F)



##### Extraction of genes with amp > 3 #####
RNA <- read.csv("RNA_stat_maxRNAover2.csv",header=T,sep=",")
RNA.over3 <- subset(RNA,ampRNA>3)
nrow(RNA.over3)
write.csv(RNA.over3, 
          file="RNA_stat_maxRNAover2_ampRNAover3.csv", 
          quote=F, row.names=F)



##### Overwriting of the file with overlap of ampover1, edgeR, and cosinor #####

# Perform "11_RNA_edgeR_cosinor.R" and "12_Ampover1_edgeR_cosinor.R" in advance

RNA <- read.csv("RNA_stat_maxRNAover2_ampRNAover3_edgeR_cosinor.csv",header=T,sep=",")
nrow(RNA)
write.csv(RNA, file="RNA_stat_maxRNAover2_ampRNAover3.csv", quote=F, row.names=F)



##### Ordering of genes according to amplitude of RNA #####
RNA <- read.csv("RNA_stat_maxRNAover2_ampRNAover3.csv",header=T,sep=",")
order <- order(RNA$ampRNA, decreasing=T)
RNA.amp.order <- RNA[order,]
write.csv(RNA.amp.order, 
          file="RNA_stat_maxRNAover2_ampRNAover3_ampRNAorder.csv",
          quote=F, row.names=F)



##### Histogram of maximum #####
RNA <- read.csv("RNA_stat.csv",header=T,sep=",")

maxRNA.mod <- rep(NA, length=nrow(RNA))
for(i in 1:nrow(RNA)){
  if(RNA$maxRNA[i] <= 8){
    maxRNA.mod[i] <- RNA$maxRNA[i]
  }else{maxRNA.mod[i] <- 8.001}
}

pdf("../figs/histogram/histogram_maxRNAmod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(maxRNA.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,8), ylim=c(0,4000),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,4000,2000),labels=F)
mtext(c("0","2","4","6","> 8"), side=1, at=seq(0,8,2),line=-0.3)
mtext(c("0","2,000","4,000"), side=2, at=seq(0,4000,2000),las=1,line=0.15)
mtext(expression(paste(log[2], "(maximum rpkm)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("RNA expression",side=3,line=0)
text(4.4,2500,"15,390",col="red")
arrows(2,0,2,4000,length=0,col="red")
arrows(2,2500,3,2500,length=0.05,col="red")
dev.off()



##### Histogram of amplitude #####
RNA <- read.csv("RNA_stat_maxRNAover2.csv",header=T,sep=",")
nrow(RNA)

ampRNA.mod <- rep(NA, length=nrow(RNA))
for(i in 1:nrow(RNA)){
  if(RNA$ampRNA[i] <= 4){
    ampRNA.mod[i] <- RNA$ampRNA[i]
  }else{ampRNA.mod[i] <- 4.001}
}

pdf("../figs/histogram/histogram_ampRNAmod.pdf",width=1.5,height=1.1)
par(ps=6)
par(mar=c(1.5,2,1,0.5))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(ampRNA.mod, las=1, tcl=-0.2, xlab="", ylab="", main="",xaxt="n", yaxt="n",
     xlim=c(0,4), ylim=c(0,1500),breaks=48,col="grey")
axis(1, tcl=-0.1, labels=F)
axis(2, tcl=-0.1, at=seq(0,1500,500),labels=F)
mtext(c("0","1","2","3","> 4"), side=1, at=seq(0,4,1),line=-0.3)
mtext(c("0","500","1,000","1,500"), side=2, at=seq(0,1500,500),las=1,line=0.15)
mtext(expression(paste(log[2], "(seasonal amplitude)")),side=1,line=0.2)
mtext("No. of genes",side=2,line=1.2)
mtext("RNA expression",side=3,line=0)
text(3.9,1000,"2,187",col="red")
arrows(3,0,3,1500,length=0,col="red")
arrows(3,1000,3.4,1000,length=0.05,col="red")
dev.off()



#####  Extraction of genes with amp > 3 for monthly mean data #####
RNA <- read.csv("RNA_stat_maxRNAover2_ampRNAover3_ampRNAorder.csv",header=T,sep=",")
RNAmean <- read.csv("RNA_monthmean.csv",header=T,sep=",")

RNA.seasonality.mean <- 
  merge(RNAmean,RNA, 
        by=c("Ahal_ID"),
        all=F,sort=F)
RNA.seasonality.mean <- RNA.seasonality.mean[,1:ncol(RNAmean)]
nrow(RNA.seasonality.mean)

write.csv(RNA.seasonality.mean, 
          file="RNA_monthmean_ampRNAover3.csv", quote=F, row.names=F)



#####  Extraction of genes with ampRNA > 3 & ampHistone > 1 #####
RNA <- read.csv("RNA_monthmean_ampRNAover3.csv",header=T,sep=",")
K4 <- read.csv("K4K27_monthmean_ampK4over1.csv",header=T,sep=",")
K27 <- read.csv("K4K27_monthmean_ampK27over1.csv",header=T,sep=",")
nrow(RNA)
nrow(K4)
nrow(K27)

RNA.K4 <- merge(RNA,K4, by=c("Ahal_ID"),all=F,sort=F)
RNA.K4 <- RNA.K4[,c(1,16:21,2:13,22:45,14:15,46:ncol(RNA.K4))]
nrow(RNA.K4)
write.csv(RNA.K4, 
          file="RNAK4K27_monthmean_ampRNAover3ampK4over1.csv", quote=F, row.names=F)

RNA.K27 <- merge(RNA,K27, by=c("Ahal_ID"),all=F,sort=F)
RNA.K27 <- RNA.K27[,c(1,16:21,2:13,22:45,14:15,46:ncol(RNA.K27))]
nrow(RNA.K27)
write.csv(RNA.K27, 
          file="RNAK4K27_monthmean_ampRNAover3ampK27over1.csv", quote=F, row.names=F)



##### Spearman's rank correlation between seasonal changes of RNA and histone modifications #####

#RNA and K4
RNAK4<-read.csv("RNAK4K27_monthmean_ampRNAover3ampK4over1.csv",header=T,sep=",")
RNA <- RNAK4[,8:19]
K4 <- RNAK4[,20:31]
RNA <- as.matrix(RNA)
K4 <- as.matrix(K4)
title <- RNAK4[,1:7]

Cor.RNA.K4 <- rep(NA,length=nrow(RNA))
cor.mat <- cbind(title, Cor.RNA.K4)
for(i in 1:nrow(RNA)){
  cor.mat[i,8] <- round(cor(RNA[i,],K4[i,], method="spearman"), digits=2)
}

sortlist <- order(as.numeric(cor.mat[,8]),decreasing=T)
cor.mat.sorted <- cor.mat[sortlist,]
write.csv(cor.mat.sorted, 
          file="RNAK4_ampRNAover3ampK4over1_cor.csv",quote=F, row.names=F)

pdf("../figs/histogram/histogram_RNAK4_ampRNAover3ampK4over1_cor.pdf",width=1.6,height=1.2)
par(ps=6)
par(mar=c(2,2,1.3,1))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(as.numeric(cor.mat.sorted[,8]), las=1, tcl=-0.2, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-1,1), ylim=c(0,150),breaks=24,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0,0),labels=F)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,150,50))
label <- c("-1.0","-0.5","0.0","0.5","1.0")
mtext(label,side=1,at=c(-1,-0.5,0,0.5,1),line=-0.35,las=1)
mtext("Correlation between",side=1,line=0.1)
mtext("seasonal dynamics of",side=1,line=0.5)
mtext("mRNA and H3K4me3",side=1,line=0.9)
mtext("No. of genes",side=2,line=0.8)
#mtext("Genes with",side=3,line=0.4)
#mtext("dual-seasonality",side=3,line=0)
arrows(0,0,0,150,length=0,col="black",lty="dashed",lwd=0.5)
dev.off()


# RNA and K27
RNAK27<-read.csv("RNAK4K27_monthmean_ampRNAover3ampK27over1.csv",header=T,sep=",")
RNA <- RNAK27[,8:19]
K27 <- RNAK27[,32:43]
RNA <- as.matrix(RNA)
K27 <- as.matrix(K27)
title <- RNAK27[,1:7]

Cor.RNA.K27 <- rep(NA,length=nrow(RNA))
cor.mat <- cbind(title, Cor.RNA.K27)
for(i in 1:nrow(RNA)){
  cor.mat[i,8] <- round(cor(RNA[i,],K27[i,], method="spearman"), digits=2)
}

sortlist <- order(as.numeric(cor.mat[,8]),decreasing=F)
cor.mat.sorted <- cor.mat[sortlist,]
write.csv(cor.mat.sorted, 
          file="RNAK27_ampRNAover3ampK27over1_cor.csv",quote=F, row.names=F)

pdf("../figs/histogram/histogram_RNAK27_ampRNAover3ampK27over1_cor.pdf",width=1.6,height=1.2)
par(ps=6)
par(mar=c(2,2,1.3,1))
par(mgp=c(0,0.3,0))
par(xpd=T)
hist(as.numeric(cor.mat.sorted[,8]), las=1, tcl=-0.2, xlab="", ylab="", main="",
     xaxt="n", yaxt="n", xlim=c(-1,1), ylim=c(0,50),breaks=24,col="grey")
axis(1, tcl=-0.1, mgp=c(0,0,0),labels=F)
axis(2, tcl=-0.1, las=1, mgp=c(0,0.15,0), at=seq(0,50,10))
label <- c("-1.0","-0.5","0.0","0.5","1.0")
mtext(label,side=1,at=c(-1,-0.5,0,0.5,1),line=-0.35,las=1)
mtext("Correlation between",side=1,line=0.1)
mtext("seasonal dynamics of",side=1,line=0.5)
mtext("mRNA and H3K27me3",side=1,line=0.9)
mtext("No. of genes",side=2,line=0.6)
#mtext("Genes with",side=3,line=0.4)
#mtext("dual-seasonality",side=3,line=0)
arrows(0,0,0,50,length=0,col="black",lty="dashed",lwd=0.5)
dev.off()



##### Mean and amplitude of RNA expression for all combinations of presence/absence of K27- and K4- enrichment and their seasonalities #####

# Load data
all <- read.csv("K4K27_data_stat.csv",header=T,sep=",")
maxK4over2<-read.csv("K4K27_data_stat_maxK4over2.csv",header=T,sep=",")
maxK27over2<-read.csv("K4K27_data_stat_maxK27over2.csv",header=T,sep=",")
maxK4K27over2<-read.csv("K4K27_data_stat_maxK4K27over2.csv",header=T,sep=",")

maxK4over2ampK4over1<-read.csv("K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
maxK27over2ampK27over1<-read.csv("K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
category9<-read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv",header=T,sep=",")

maxK4over2ampK4over1ampK27under1<-read.csv("K4K27_data_stat_maxK4over2_ampK4onlyover1.csv",header=T,sep=",")
maxK27over2ampK4under1ampK27over1<-read.csv("K4K27_data_stat_maxK27over2_ampK27onlyover1.csv",header=T,sep=",")

# Extraction of genes classified into each category

#category         1    2    3    4    5    6    7    8    9
#K27enriched           yes  yes            yes  yes  yes  yes
#K27seasonality             yes                 yes       yes
#K4enriched                      yes  yes  yes  yes  yes  yes
#K4seasonality                        yes            yes  yes

temp <- merge(maxK27over2,maxK4over2, by="Ahal_ID",all=T,sort=F)
category2.3 <- temp[is.na(temp$K4_rep1_start.y),1:ncol(maxK27over2)]
temp <- merge(category2.3,maxK27over2ampK27over1, by="Ahal_ID",all=T,sort=F)
category2 <- temp[is.na(temp$K4_rep1_start),1:ncol(category2.3)]
temp <- merge(maxK27over2ampK27over1,maxK4over2, by="Ahal_ID",all=T,sort=F)
category3 <- temp[is.na(temp$K4_rep1_start.y),1:ncol(maxK27over2ampK27over1)]

temp <- merge(maxK4over2,maxK27over2, by="Ahal_ID",all=T,sort=F)
category4.5 <- temp[is.na(temp$K4_rep1_start.y),1:ncol(maxK4over2)]
temp <- merge(category4.5,maxK4over2ampK4over1, by="Ahal_ID",all=T,sort=F)
category4 <- temp[is.na(temp$K4_rep1_start),1:ncol(category4.5)]
temp <- merge(maxK4over2ampK4over1,maxK27over2, by="Ahal_ID",all=T,sort=F)
category5 <- temp[is.na(temp$K4_rep1_start.y),1:ncol(maxK4over2ampK4over1)]

category8.9 <- merge(maxK4over2ampK4over1,maxK27over2, by="Ahal_ID",all=F,sort=F)
temp <- merge(category8.9,category9, by="Ahal_ID",all=T,sort=F)
category8 <- temp[is.na(temp$K4_rep1_start),1:ncol(category9)]

category7.9 <- merge(maxK27over2ampK27over1,maxK4over2, by="Ahal_ID",all=F,sort=F)
temp <- merge(category7.9,category9, by="Ahal_ID",all=T,sort=F)
category7 <- temp[is.na(temp$K4_rep1_start),1:ncol(category9)]

temp <- merge(maxK4K27over2,category8.9, by="Ahal_ID",all=T,sort=F)
category6.7 <- temp[is.na(temp$K4_rep1_start.x),1:ncol(maxK4K27over2)]
temp <- merge(category6.7,category7, by="Ahal_ID",all=T,sort=F)
category6 <- temp[is.na(temp$K4_rep1_start.x),1:ncol(category6.7)]

temp <- merge(all,maxK27over2, by="Ahal_ID",all=T,sort=F)
category1.4.5 <- temp[is.na(temp$K4_rep1_start.y),1:ncol(all)]
temp <- merge(category1.4.5,category4.5, by="Ahal_ID",all=T,sort=F)
category1 <- temp[is.na(temp$K4_rep1_start.x.y),1:ncol(category1.4.5)]

print(c(nrow(category1),nrow(category2),nrow(category3),
    nrow(category4),nrow(category5),nrow(category6),
    nrow(category7),nrow(category8),nrow(category9)))

sum(nrow(category1),nrow(category2),nrow(category3),
    nrow(category4),nrow(category5),nrow(category6),
    nrow(category7),nrow(category8),nrow(category9))
    

colnames <- vector(length=ncol(category1))
list <- strsplit(colnames(category1), "\\.x\\.x")
for(i in 1:ncol(category1)){
  colnames[i] <- list[[i]][1]
}
colnames(category1) <- colnames


category.list <-list(category1,category2,category3,category4,category5,category6,category7,category8,category9)

for(i in c(2,3,4,5,7,8)){
	colnames <- vector(length=ncol(category.list[[i]]))
	list <- strsplit(colnames(category.list[[i]]), "\\.x")
	for(j in 1:ncol(category.list[[i]])){
  	colnames[j] <- list[[j]][1]
	}
	colnames(category.list[[i]]) <- colnames
}

for(i in 1:9){
	write.csv(category.list[[i]], 
    	      file=paste0("K4K27_data_stat_category",i,".csv"), quote=F, row.names=F)
}


# Combine each category with RNA.stat
RNA.stat<-read.csv("RNA_stat.csv",header=T,sep=",")
category1<-read.csv("K4K27_data_stat_category1.csv",header=T,sep=",")
category2<-read.csv("K4K27_data_stat_category2.csv",header=T,sep=",")
category3<-read.csv("K4K27_data_stat_category3.csv",header=T,sep=",")
category4<-read.csv("K4K27_data_stat_category4.csv",header=T,sep=",")
category5<-read.csv("K4K27_data_stat_category5.csv",header=T,sep=",")
category6<-read.csv("K4K27_data_stat_category6.csv",header=T,sep=",")
category7<-read.csv("K4K27_data_stat_category7.csv",header=T,sep=",")
category8<-read.csv("K4K27_data_stat_category8.csv",header=T,sep=",")
category9<-read.csv("K4K27_data_stat_category9.csv",header=T,sep=",")

temp <- merge(RNA.stat,category1,by="Ahal_ID",all=F,sort=F)
RNA.stat.category1 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category2,by="Ahal_ID",all=F,sort=F)
RNA.stat.category2 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category3,by="Ahal_ID",all=F,sort=F)
RNA.stat.category3 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category4,by="Ahal_ID",all=F,sort=F)
RNA.stat.category4 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category5,by="Ahal_ID",all=F,sort=F)
RNA.stat.category5 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category6,by="Ahal_ID",all=F,sort=F)
RNA.stat.category6 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category7,by="Ahal_ID",all=F,sort=F)
RNA.stat.category7 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category8,by="Ahal_ID",all=F,sort=F)
RNA.stat.category8 <- temp[,1:ncol(RNA.stat)]

temp <- merge(RNA.stat,category9,by="Ahal_ID",all=F,sort=F)
RNA.stat.category9 <- temp[,1:ncol(RNA.stat)]


# Wilcoxon rank sum test
# K27 seasonality
RNAmean_2_3 <- wilcox.exact(x=RNA.stat.category2$meanRNA, y=RNA.stat.category3$meanRNA, paired=F)
RNAmean_6_7 <- wilcox.exact(x=RNA.stat.category6$meanRNA, y=RNA.stat.category7$meanRNA, paired=F)
RNAmean_8_9 <- wilcox.exact(x=RNA.stat.category8$meanRNA, y=RNA.stat.category9$meanRNA, paired=F)

RNAamp_2_3 <-wilcox.exact(x=RNA.stat.category2$ampRNA, y=RNA.stat.category3$ampRNA, paired=F)
RNAamp_6_7 <- wilcox.exact(x=RNA.stat.category6$ampRNA, y=RNA.stat.category7$ampRNA, paired=F)
RNAamp_8_9 <- wilcox.exact(x=RNA.stat.category8$ampRNA, y=RNA.stat.category9$ampRNA, paired=F)

# K4 seasonality
RNAmean_4_5 <- wilcox.exact(x=RNA.stat.category4$meanRNA, y=RNA.stat.category5$meanRNA, paired=F)
RNAmean_6_8 <- wilcox.exact(x=RNA.stat.category6$meanRNA, y=RNA.stat.category8$meanRNA, paired=F)
RNAmean_7_9 <- wilcox.exact(x=RNA.stat.category7$meanRNA, y=RNA.stat.category9$meanRNA, paired=F)

RNAamp_4_5 <- wilcox.exact(x=RNA.stat.category4$ampRNA, y=RNA.stat.category5$ampRNA, paired=F)
RNAamp_6_8 <- wilcox.exact(x=RNA.stat.category6$ampRNA, y=RNA.stat.category8$ampRNA, paired=F)
RNAamp_7_9 <- wilcox.exact(x=RNA.stat.category7$ampRNA, y=RNA.stat.category9$ampRNA, paired=F)

adp.mean <- ng.BHFDR(c(RNAmean_2_3$p.value,RNAmean_6_7$p.value,RNAmean_8_9$p.value,
				  RNAmean_4_5$p.value,RNAmean_6_8$p.value,RNAmean_7_9$p.value))
adp.mean

adp.amp <- ng.BHFDR(c(RNAamp_2_3$p.value,RNAamp_6_7$p.value,RNAamp_8_9$p.value,
				  RNAamp_4_5$p.value,RNAamp_6_8$p.value,RNAamp_7_9$p.value))
adp.amp


# Brunner-Munzel test

# K27 seasonality
RNAmean_2_3 <- brunner.munzel.test(x=RNA.stat.category2$meanRNA, y=RNA.stat.category3$meanRNA,alternative="less")
RNAmean_6_7 <- brunner.munzel.test(x=RNA.stat.category6$meanRNA, y=RNA.stat.category7$meanRNA,alternative="less")
RNAmean_8_9 <- brunner.munzel.test(x=RNA.stat.category8$meanRNA, y=RNA.stat.category9$meanRNA,alternative="less")

RNAamp_2_3 <-brunner.munzel.test(x=RNA.stat.category2$ampRNA, y=RNA.stat.category3$ampRNA,alternative="less")
RNAamp_6_7 <- brunner.munzel.test(x=RNA.stat.category6$ampRNA, y=RNA.stat.category7$ampRNA,alternative="less")
RNAamp_8_9 <- brunner.munzel.test(x=RNA.stat.category8$ampRNA, y=RNA.stat.category9$ampRNA,alternative="less")

# K4 seasonality
RNAmean_4_5 <- brunner.munzel.test(x=RNA.stat.category4$meanRNA, y=RNA.stat.category5$meanRNA,alternative="greater")
RNAmean_6_8 <- brunner.munzel.test(x=RNA.stat.category6$meanRNA, y=RNA.stat.category8$meanRNA,alternative="greater")
RNAmean_7_9 <- brunner.munzel.test(x=RNA.stat.category7$meanRNA, y=RNA.stat.category9$meanRNA,alternative="greater")

RNAamp_4_5 <- brunner.munzel.test(x=RNA.stat.category4$ampRNA, y=RNA.stat.category5$ampRNA,alternative="less")
RNAamp_6_8 <- brunner.munzel.test(x=RNA.stat.category6$ampRNA, y=RNA.stat.category8$ampRNA,alternative="less")
RNAamp_7_9 <- brunner.munzel.test(x=RNA.stat.category7$ampRNA, y=RNA.stat.category9$ampRNA,alternative="less")


adp.mean.K27 <- ng.BHFDR(c(RNAmean_2_3$p.value,RNAmean_6_7$p.value,RNAmean_8_9$p.value))
adp.mean.K27

adp.mean.K4 <- ng.BHFDR(c(RNAmean_4_5$p.value,RNAmean_6_8$p.value,RNAmean_7_9$p.value))
adp.mean.K4

adp.amp.K27 <- ng.BHFDR(c(RNAamp_2_3$p.value,RNAamp_6_7$p.value,RNAamp_8_9$p.value))
adp.amp.K27

adp.amp.K4 <- ng.BHFDR(c(RNAamp_4_5$p.value,RNAamp_6_8$p.value,RNAamp_7_9$p.value))
adp.amp.K4
				  

# Steel-Dwass test
data <- 
  c(RNA.stat.category1$meanRNA, RNA.stat.category2$meanRNA, RNA.stat.category3$meanRNA,
    RNA.stat.category4$meanRNA, RNA.stat.category5$meanRNA, RNA.stat.category6$meanRNA, RNA.stat.category7$meanRNA, RNA.stat.category8$meanRNA, RNA.stat.category9$meanRNA)
data2 <- 
  c(RNA.stat.category1$ampRNA, RNA.stat.category2$ampRNA, RNA.stat.category3$ampRNA,
    RNA.stat.category4$ampRNA, RNA.stat.category5$ampRNA, RNA.stat.category6$ampRNA, RNA.stat.category7$ampRNA, RNA.stat.category8$ampRNA, RNA.stat.category9$ampRNA)
group <- rep(1:9, c(nrow(category1),nrow(category2),nrow(category3),nrow(category4),nrow(category5),nrow(category6),nrow(category7),nrow(category8),nrow(category9)))

source("http://aoki2.si.gunma-u.ac.jp/R/src/Steel-Dwass.R", encoding="euc-jp")
Steel.Dwass(data,group)
Steel.Dwass(data2,group)



##### Visualization by boxplot #####
pdf("../figs/boxplot/boxplot_meanRNA_allcategories.pdf",
    width=2,height=1.5)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(3,4.6,1,0.4))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,9.5),ylim=c(0,8))
axis(1,at=1:9,labels=F)
axis(2,at=seq(0,8,2),las=1, labels=F)
b <-
boxplot(RNA.stat.category1$meanRNA, RNA.stat.category2$meanRNA, RNA.stat.category3$meanRNA,
        RNA.stat.category4$meanRNA, RNA.stat.category5$meanRNA, RNA.stat.category6$meanRNA,
        RNA.stat.category7$meanRNA, RNA.stat.category8$meanRNA, RNA.stat.category9$meanRNA,
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.14,-0.32,9.5,-0.32,length=0,lwd=1)
arrows(0.14,-0.32,0.14,0,length=0,lwd=1)

text(1,b$stats[5,1]+0.6,"e",cex=0.8)
text(2,b$stats[5,2]+0.6,"f",cex=0.8)
text(3,b$stats[5,3]+0.6,"f",cex=0.8)
text(4,b$stats[5,4]+0.6,"a",cex=0.8)
text(5,b$stats[5,5]+0.6,"b",cex=0.8)
text(6,b$stats[5,6]+0.6,"d",cex=0.8)
text(7,b$stats[5,7]+0.6,"bc",cex=0.8)
text(8,b$stats[5,8]+0.6,"c",cex=0.8)
text(9,b$stats[5,9]+0.6,"bc",cex=0.8)

points(c(2:3,6:9),rep(-1.5,6),pch=3,cex=0.4)
points(c(3,7,9),rep(-2.5,3),pch=3,cex=0.4)
points(4:9,rep(-3.5,6),pch=3,cex=0.4)
points(c(5,8,9),rep(-4.5,3),pch=3,cex=0.4)
text(-3.92,-1.5,"H3K27me3 enrichment")
text(-3.9,-2.5,"H3K27me3 seasonality")
text(-3.72,-3.5,"H3K4me3 enrichment")
text(-3.7,-4.5,"H3K4me3 seasonality")
mtext(c("0","2","4","6","8"), side=2, line=0.15, at=seq(0,8,2), las=1)
mtext(expression(paste(log[2],"(rpkm)")), side=2, line=0.3)
mtext("Mean expression", side=3, line=0)
dev.off()


pdf("../figs/boxplot/boxplot_ampRNA_allcategories.pdf",
    width=2,height=1.5)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(3,4.6,1,0.4))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,9.5),ylim=c(0,8))
axis(1,at=1:9,labels=F)
axis(2,at=seq(0,8,2),las=1, labels=F)
b <-
boxplot(RNA.stat.category1$ampRNA, RNA.stat.category2$ampRNA, RNA.stat.category3$ampRNA,
        RNA.stat.category4$ampRNA, RNA.stat.category5$ampRNA, RNA.stat.category6$ampRNA,
        RNA.stat.category7$ampRNA, RNA.stat.category8$ampRNA, RNA.stat.category9$ampRNA,
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.14,-0.32,9.5,-0.32,length=0,lwd=1)
arrows(0.14,-0.32,0.14,0,length=0,lwd=1)

text(1,b$stats[5,1]+0.6,"e",cex=0.8)
text(2,b$stats[5,2]+0.6,"f",cex=0.8)
text(3,b$stats[5,3]+0.6,"f",cex=0.8)
text(4,b$stats[5,4]+0.6,"b",cex=0.8)
text(5,b$stats[5,5]+0.6,"b",cex=0.8)
text(6,b$stats[5,6]+0.6,"d",cex=0.8)
text(7,b$stats[5,7]+0.6,"c",cex=0.8)
text(8,b$stats[5,8]+0.6,"bc",cex=0.8)
text(9,b$stats[5,9]+0.6,"a",cex=0.8)

points(c(2:3,6:9),rep(-1.5,6),pch=3,cex=0.4)
points(c(3,7,9),rep(-2.5,3),pch=3,cex=0.4)
points(4:9,rep(-3.5,6),pch=3,cex=0.4)
points(c(5,8,9),rep(-4.5,3),pch=3,cex=0.4)
text(-3.92,-1.5,"H3K27me3 enrichment")
text(-3.9,-2.5,"H3K27me3 seasonality")
text(-3.72,-3.5,"H3K4me3 enrichment")
text(-3.7,-4.5,"H3K4me3 seasonality")
mtext(c("0","2","4","6","8"), side=2, line=0.15, at=seq(0,8,2), las=1)
mtext(expression(paste(log[2],"(rpkm)")), side=2, line=0.3)
mtext("Seasonal amplitude of", side=3, line=0.3)
mtext("expression", side=3, line=-0.1)
dev.off()


# H3K27me3seasonality 
pdf("../figs/boxplot/boxplot_meanRNA_K27seasonality.pdf",
    width=1.8,height=1.7)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(3,4.4,1,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,6.5),ylim=c(0,8))
axis(1,at=1:6,labels=F)
axis(2,at=seq(0,8,2),las=1, labels=F)
b <-
boxplot(RNA.stat.category2$meanRNA, RNA.stat.category3$meanRNA, 
        RNA.stat.category6$meanRNA, RNA.stat.category7$meanRNA,
        RNA.stat.category8$meanRNA, RNA.stat.category9$meanRNA,
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.26,-0.32,6.5,-0.32,length=0,lwd=1)
arrows(0.26,-0.32,0.26,0,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.3,1,1,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.3,2,1,length=0,lwd=0.5)
arrows(1,1,2,1,length=0,lwd=0.5)
arrows(3,b$stats[5,3]+0.3,3,6.9,length=0,lwd=0.5)
arrows(4,b$stats[5,4]+0.3,4,6.9,length=0,lwd=0.5)
arrows(3,6.9,4,6.9,length=0,lwd=0.5)
arrows(5,b$stats[5,5]+0.3,5,5.8,length=0,lwd=0.5)
arrows(6,b$stats[5,6]+0.3,6,5.8,length=0,lwd=0.5)
arrows(5,5.8,6,5.8,length=0,lwd=0.5)
text(1.5,1.4,"n.s.",cex=0.8)
text(2.5,7.4,expression(paste(italic(P)," = ","9.6 × ",10^-12)),cex=0.8)
text(5.2,6.9,expression(paste(italic(P)," =",)),cex=0.8)
text(6,6.3,expression(paste("2.3 × ",10^-3)),cex=0.8)

points(1:6,rep(-1.2,6),pch=3,cex=0.4)
points(c(2,4,6),rep(-2.1,3),pch=3,cex=0.4)
points(3:6,rep(-3,4),pch=3,cex=0.4)
points(5:6,rep(-3.9,2),pch=3,cex=0.4)
text(-3.12,-1.2,"H3K27me3 enrichment")
text(-3.1,-2.1,"H3K27me3 seasonality")
text(-2.92,-3,"H3K4me3 enrichment")
text(-2.9,-3.9,"H3K4me3 seasonality")
mtext(c("0","2","4","6","8"), side=2, line=0.15, at=seq(0,8,2), las=1)
mtext(expression(paste(log[2],"(rpkm)")), side=2, line=0.3)
mtext("Mean expression", side=3, line=0.1)
dev.off()


pdf("../figs/boxplot/boxplot_ampRNA_K27seasonality.pdf",
    width=1.8,height=1.7)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(3,4.4,1,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,6.5),ylim=c(0,8))
axis(1,at=1:6,labels=F)
axis(2,at=seq(0,8,2),las=1, labels=F)
b <-
boxplot(RNA.stat.category2$ampRNA, RNA.stat.category3$ampRNA, 
        RNA.stat.category6$ampRNA, RNA.stat.category7$ampRNA,
        RNA.stat.category8$ampRNA, RNA.stat.category9$ampRNA,
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.26,-0.32,6.5,-0.32,length=0,lwd=1)
arrows(0.26,-0.32,0.26,0,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.3,1,1.8,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.3,2,1.8,length=0,lwd=0.5)
arrows(1,1.8,2,1.8,length=0,lwd=0.5)
arrows(3,b$stats[5,3]+0.3,3,5.8,length=0,lwd=0.5)
arrows(4,b$stats[5,4]+0.3,4,5.8,length=0,lwd=0.5)
arrows(3,5.8,4,5.8,length=0,lwd=0.5)
arrows(5,b$stats[5,5]+0.3,5,7.9,length=0,lwd=0.5)
arrows(6,b$stats[5,6]+0.3,6,7.9,length=0,lwd=0.5)
arrows(5,7.9,6,7.9,length=0,lwd=0.5)
text(1.5,2.2,"n.s.",cex=0.8)
text(2.5,6.3,expression(paste(italic(P)," = ","4.3 × ",10^-15)),cex=0.8)
text(3.2,7.7,expression(paste(italic(P)," = ","7.9 × ",10^-6)),cex=0.8)

#text(5.2,6.8,expression(paste(italic(P)," <",)),cex=0.8)
#text(5.5,6.2,"0.001",cex=0.8)

points(1:6,rep(-1.2,6),pch=3,cex=0.4)
points(c(2,4,6),rep(-2.1,3),pch=3,cex=0.4)
points(3:6,rep(-3,4),pch=3,cex=0.4)
points(5:6,rep(-3.9,2),pch=3,cex=0.4)
text(-3.12,-1.2,"H3K27me3 enrichment")
text(-3.1,-2.1,"H3K27me3 seasonality")
text(-2.92,-3,"H3K4me3 enrichment")
text(-2.9,-3.9,"H3K4me3 seasonality")
mtext(c("0","2","4","6","8"), side=2, line=0.15, at=seq(0,8,2), las=1)
mtext(expression(paste(log[2],"(rpkm)")), side=2, line=0.3)
mtext("Seasonal amplitude of", side=3, line=0.4)
mtext("expression", side=3, line=0)
dev.off()


# H3K4me3seasonality 
pdf("../figs/boxplot/boxplot_meanRNA_K4seasonality.pdf",
    width=1.8,height=1.7)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(3,4.4,1,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,6.5),ylim=c(0,8))
axis(1,at=1:6,labels=F)
axis(2,at=seq(0,8,2),las=1, labels=F)
b <- 
boxplot(RNA.stat.category4$meanRNA, RNA.stat.category5$meanRNA, 
        RNA.stat.category6$meanRNA, RNA.stat.category8$meanRNA,
        RNA.stat.category7$meanRNA, RNA.stat.category9$meanRNA,
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.26,-0.32,6.5,-0.32,length=0,lwd=1)
arrows(0.26,-0.32,0.26,1,length=0,lwd=1)
        
arrows(1,b$stats[5,1]+0.3,1,8.2,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.3,2,8.2,length=0,lwd=0.5)
arrows(1,8.2,2,8.2,length=0,lwd=0.5)
arrows(3,b$stats[5,3]+0.3,3,5.8,length=0,lwd=0.5)
arrows(4,b$stats[5,4]+0.3,4,5.8,length=0,lwd=0.5)
arrows(3,5.8,4,5.8,length=0,lwd=0.5)
arrows(5,b$stats[5,5]+0.3,5,6.8,length=0,lwd=0.5)
arrows(6,b$stats[5,6]+0.3,6,6.8,length=0,lwd=0.5)
arrows(5,6.8,6,6.8,length=0,lwd=0.5)
text(3.9,8.2,expression(paste(italic(P)," < ","1.0 × ",10^-15)),cex=0.8)
text(3.5,6.2,"n.s.",cex=0.8)
text(5.5,7.2,"n.s.",cex=0.8)

points(3:6,rep(-1.2,4),pch=3,cex=0.4)
points(5:6,rep(-2.1,2),pch=3,cex=0.4)
points(1:6,rep(-3,6),pch=3,cex=0.4)
points(c(2,4,6),rep(-3.9,3),pch=3,cex=0.4)
text(-3.12,-1.2,"H3K27me3 enrichment")
text(-3.1,-2.1,"H3K27me3 seasonality")
text(-2.92,-3,"H3K4me3 enrichment")
text(-2.9,-3.9,"H3K4me3 seasonality")
mtext(c("0","2","4","6","8"), side=2, line=0.15, at=seq(0,8,2), las=1)
mtext(expression(paste(log[2],"(rpkm)")), side=2, line=0.3)
mtext("Mean expression", side=3, line=0.1)
dev.off()


pdf("../figs/boxplot/boxplot_ampRNA_K4seasonality.pdf",
    width=1.8,height=1.7)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(3,4.4,1,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,6.5),ylim=c(0,8))
axis(1,at=1:6,labels=F)
axis(2,at=seq(0,8,2),las=1, labels=F)
b <-
boxplot(RNA.stat.category4$ampRNA, RNA.stat.category5$ampRNA, 
        RNA.stat.category6$ampRNA, RNA.stat.category8$ampRNA,
        RNA.stat.category7$ampRNA, RNA.stat.category9$ampRNA,
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.26,-0.32,6.5,-0.32,length=0,lwd=1)
arrows(0.26,-0.32,0.26,1,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.3,1,6.4,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.3,2,6.4,length=0,lwd=0.5)
arrows(1,6.4,2,6.4,length=0,lwd=0.5)
arrows(3,b$stats[5,3]+0.3,3,6.8,length=0,lwd=0.5)
arrows(4,b$stats[5,4]+0.3,4,6.8,length=0,lwd=0.5)
arrows(3,6.8,4,6.8,length=0,lwd=0.5)
arrows(5,b$stats[5,5]+0.3,5,7.9,length=0,lwd=0.5)
arrows(6,b$stats[5,6]+0.3,6,7.9,length=0,lwd=0.5)
arrows(5,7.9,6,7.9,length=0,lwd=0.5)
text(1.5,6.8,"n.s.",cex=0.8)
text(2.7,7.7,expression(paste(italic(P)," < ")),cex=0.8)
text(3.4,7.2,expression(paste("1.0 × ",10^-15)),cex=0.8)
text(5.6,8.3,expression(paste(italic(P)," = ","1.1 × ",10^-11)),cex=0.8)

points(3:6,rep(-1.2,4),pch=3,cex=0.4)
points(5:6,rep(-2.1,2),pch=3,cex=0.4)
points(1:6,rep(-3,6),pch=3,cex=0.4)
points(c(2,4,6),rep(-3.9,3),pch=3,cex=0.4)
text(-3.12,-1.2,"H3K27me3 enrichment")
text(-3.1,-2.1,"H3K27me3 seasonality")
text(-2.92,-3,"H3K4me3 enrichment")
text(-2.9,-3.9,"H3K4me3 seasonality")
mtext(c("0","2","4","6","8"), side=2, line=0.15, at=seq(0,8,2), las=1)
mtext(expression(paste(log[2],"(rpkm)")), side=2, line=0.3)
mtext("Seasonal amplitude of", side=3, line=0.4)
mtext("expression", side=3, line=0)
dev.off()
