
library(fields)
source("functions/GOanalysis_functions.R")



##### Correlation with climate #####
sma.period <- c(7,30)

# H3K27me3
K27 <- read.csv("data/K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
K27.rep1.rep2 <- K27[,c(26:37,56:67)]
head(K27.rep1.rep2)
title <- K27[,1:7]

climate <- read.csv("data/SMA_climate.csv",header=T,sep=",")

sum.table <- matrix(nrow=nrow(K27),ncol=8)
for (i in 1:nrow(K27)){
  for (j in 1:8){
  y <- as.numeric(K27.rep1.rep2[i,])
  x <- rep(climate[,j+1],2)
  sum.table[i,j] <- cor(x,y,method="spearman")
  }
}

colnames(sum.table) <- paste0("cor.K27.", 
                             c(paste0("SMA",sma.period,c("_temp")),
                               paste0("SMA",sma.period,c("_ppt")),
                               paste0("SMA",sma.period,c("_day")),
                               paste0("SMA",sma.period,c("_sun")))
)

sum.table2 <- cbind(title,sum.table)
write.csv(sum.table2,"data/SMA_K27_cor_result.csv",quote=F,row.names=F)


# H3K4me3
K4 <- read.csv("data/K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
K4.rep1.rep2 <- K4[,c(11:22,41:52)]
head(K4.rep1.rep2)
title <- K4[,1:7]

climate <- read.csv("data/SMA_climate.csv",header=T,sep=",")

sum.table <- matrix(nrow=nrow(K4),ncol=8)
for (i in 1:nrow(K4)){
  for (j in 1:8){
    y <- as.numeric(K4.rep1.rep2[i,])
    x <- rep(climate[,j+1],2)
    sum.table[i,j] <- cor(x,y,method="spearman")
  }
}

colnames(sum.table) <- paste0("cor.K4.", 
                              c(paste0("SMA",sma.period,c("_temp")),
                                paste0("SMA",sma.period,c("_ppt")),
                                paste0("SMA",sma.period,c("_day")),
                                paste0("SMA",sma.period,c("_sun")))
)

sum.table2 <- cbind(title,sum.table)
write.csv(sum.table2,"data/SMA_K4_cor_result.csv",quote=F,row.names=F)



##### P value of correlation #####

# H3K27me3
K27 <- read.csv("data/K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
K27.rep1.rep2 <- K27[,c(26:37,56:67)]
head(K27.rep1.rep2)
title <- K27[,1:7]

climate <- read.csv("data/SMA_climate.csv",header=T,sep=",")


sum.table <- matrix(nrow=nrow(K27),ncol=8)
for (i in 1:nrow(K27)){
  for (j in 1:8){
    y <- as.numeric(K27.rep1.rep2[i,])
    x <- rep(climate[,j+1],2)
    sum.table[i,j] <- cor.test(x,y,method="spearman")$p.value
  }
}

colnames(sum.table) <- paste0("cor.K27.", 
                              c(paste0("SMA",sma.period,c("_temp")),
                                paste0("SMA",sma.period,c("_ppt")),
                                paste0("SMA",sma.period,c("_day")),
                                paste0("SMA",sma.period,c("_sun")))
)

sum.table2 <- cbind(title,sum.table)

write.csv(sum.table2,"data/SMA_K27_cor_pval_result.csv",quote=F,row.names=F)


# H3K4me3
K4 <- read.csv("data/K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
K4.rep1.rep2 <- K4[,c(11:22,41:52)]
head(K4.rep1.rep2)
title <- K4[,1:7]

climate <- read.csv("data/SMA_climate.csv",header=T,sep=",")


sum.table <- matrix(nrow=nrow(K4),ncol=8)
for (i in 1:nrow(K4)){
  for (j in 1:8){
    y <- as.numeric(K4.rep1.rep2[i,])
    x <- rep(climate[,j+1],2)
    sum.table[i,j] <- cor.test(x,y,method="spearman")$p.value
  }
}

colnames(sum.table) <- paste0("cor.K4.", 
                              c(paste0("SMA",sma.period,c("_temp")),
                                paste0("SMA",sma.period,c("_ppt")),
                                paste0("SMA",sma.period,c("_day")),
                                paste0("SMA",sma.period,c("_sun")))
)

sum.table2 <- cbind(title,sum.table)
write.csv(sum.table2,"data/SMA_K4_cor_pval_result.csv",quote=F,row.names=F)



##### Counting environmental terms #####

# K27
sma_K27 <- read.csv("data/SMA_K27_cor_result.csv",header=T,sep=",")
term <- c("temp","ppt","day","sun")

# SMA7
sma_K27_SMA7 <- sma_K27[,c(8,10,12,14)]

res_term_K27_SMA7 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  res_term_K27_SMA7[i] <- term[which.max(abs(sma_K27_SMA7[i,]))]
}

table(res_term_K27_SMA7)
plot(table(res_term_K27_SMA7))


# SMA30
sma_K27_SMA30 <- sma_K27[,c(9,11,13,15)]

res_term_K27_SMA30 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  res_term_K27_SMA30[i] <- term[which.max(abs(sma_K27_SMA30[i,]))]
}

table(res_term_K27_SMA30)
plot(table(res_term_K27_SMA30))


# K4
sma_K4 <- read.csv("data/SMA_K4_cor_result.csv",header=T,sep=",")
term <- c("temp","ppt","day","sun")

# SMA7
sma_K4_SMA7 <- sma_K4[,c(8,10,12,14)]

res_term_K4_SMA7 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  res_term_K4_SMA7[i] <- term[which.max(abs(sma_K4_SMA7[i,]))]
}

table(res_term_K4_SMA7)
plot(table(res_term_K4_SMA7))


# SMA30
sma_K4_SMA30 <- sma_K4[,c(9,11,13,15)]

res_term_K4_SMA30 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  res_term_K4_SMA30[i] <- term[which.max(abs(sma_K4_SMA30[i,]))]
}

table(res_term_K4_SMA30)
plot(table(res_term_K4_SMA30))


pdf("figs/climate_correlation/SMA7_Histone_cor_result.pdf",height=1.2,width=2)
set.panel(1,2)
par(mar=c(2.2,2,1,0.5))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)

names=c("Temperature","Precipitation","Day length","Sunlight hours")

barplot=barplot(table(res_term_K27_SMA7)[c("temp","ppt","day","sun")],ylim=c(0,2000),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,2000,1000), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-730,-730,-655,-820),labels=names,srt=45)
mtext(c("0","1,000","2,000"),at=seq(0,2000,1000),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K27me3",side=3,line=0.1)

barplot=barplot(table(res_term_K4_SMA7)[c("temp","ppt","day","sun")],ylim=c(0,2000),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,2000,1000), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-730,-730,-655,-820),labels=names,srt=45)
mtext(c("0","1,000","2,000"),at=seq(0,2000,1000),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K4me3",side=3,line=0.1)

dev.off()


pdf("figs/climate_correlation/SMA30_Histone_cor_result.pdf",height=1.2,width=2)
set.panel(1,2)
par(mar=c(2.2,2,1,0.5))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)

names=c("Temperature","Precipitation","Day length","Sunlight hours")

barplot=barplot(table(res_term_K27_SMA30)[c("temp","ppt","day","sun")],ylim=c(0,2000),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,2000,1000), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-730,-730,-655,-820),labels=names,srt=45)
mtext(c("0","1,000","2,000"),at=seq(0,2000,1000),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K27me3",side=3,line=0.1)

barplot=barplot(table(res_term_K4_SMA30)[c("temp","ppt","day","sun")],ylim=c(0,2000),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,2000,1000), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-730,-730,-655,-820),labels=names,srt=45)
mtext(c("0","1,000","2,000"),at=seq(0,2000,1000),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K4me3",side=3,line=0.1)

dev.off()



##### Counting environmental terms p-value < 0.001 #####

# K27
sma_K27 <- read.csv("data/SMA_K27_cor_result.csv",header=T,sep=",")
sma_K27_pval <- read.csv("data/SMA_K27_cor_pval_result.csv",header=T,sep=",")
term <- c("temp","ppt","day","sun")

# SMA7
sma_K27_SMA7 <- sma_K27[,c(8,10,12,14)]
sma_K27_pval_SMA7 <- sma_K27_pval[,c(8,10,12,14)]

res_term_K27_SMA7 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  res_term_K27_SMA7[i] <- which.max(abs(sma_K27_SMA7[i,]))
}

pval_SMA7 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  pval_SMA7[i] <- sma_K27_pval_SMA7[i,res_term_K27_SMA7[i]]
}

adp_SMA7 <- ng.BHFDR(pval_SMA7)

res_term2_K27_SMA7 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  if(adp_SMA7[i] < 0.001){
    res_term2_K27_SMA7[i] <- term[res_term_K27_SMA7[i]]
  }
}

table(res_term2_K27_SMA7)
plot(table(res_term2_K27_SMA7))
sum(table(res_term2_K27_SMA7))/nrow(sma_K27)


# SMA30
sma_K27_SMA30 <- sma_K27[,c(9,11,13,15)]
sma_K27_pval_SMA30 <- sma_K27_pval[,c(9,11,13,15)]

res_term_K27_SMA30 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  res_term_K27_SMA30[i] <- which.max(abs(sma_K27_SMA30[i,]))
}

pval_SMA30 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  pval_SMA30[i] <- sma_K27_pval_SMA30[i,res_term_K27_SMA30[i]]
}

adp_SMA30 <- ng.BHFDR(pval_SMA30)

res_term2_K27_SMA30 <- rep(NA,length=nrow(sma_K27))
for(i in 1:nrow(sma_K27)){
  if(adp_SMA30[i] < 0.001){
    res_term2_K27_SMA30[i] <- term[res_term_K27_SMA30[i]]
  }
}

table(res_term2_K27_SMA30)
plot(table(res_term2_K27_SMA30))
sum(table(res_term2_K27_SMA30))/nrow(sma_K27)


# K4
sma_K4 <- read.csv("data/SMA_K4_cor_result.csv",header=T,sep=",")
sma_K4_pval <- read.csv("data/SMA_K4_cor_pval_result.csv",header=T,sep=",")
term <- c("temp","ppt","day","sun")


# SMA7
sma_K4_SMA7 <- sma_K4[,c(8,10,12,14)]
sma_K4_pval_SMA7 <- sma_K4_pval[,c(8,10,12,14)]

res_term_K4_SMA7 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  res_term_K4_SMA7[i] <- which.max(abs(sma_K4_SMA7[i,]))
}

pval_SMA7 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  pval_SMA7[i] <- sma_K4_pval_SMA7[i,res_term_K4_SMA7[i]]
}

adp_SMA7 <- ng.BHFDR(pval_SMA7)

res_term2_K4_SMA7 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  if(adp_SMA7[i] < 0.001){
    res_term2_K4_SMA7[i] <- term[res_term_K4_SMA7[i]]
  }
}

table(res_term2_K4_SMA7)
plot(table(res_term2_K4_SMA7))
sum(table(res_term2_K4_SMA7))/nrow(sma_K4)


# SMA30
sma_K4_SMA30 <- sma_K4[,c(9,11,13,15)]
sma_K4_pval_SMA30 <- sma_K4_pval[,c(9,11,13,15)]

res_term_K4_SMA30 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  res_term_K4_SMA30[i] <- which.max(abs(sma_K4_SMA30[i,]))
}

pval_SMA30 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  pval_SMA30[i] <- sma_K4_pval_SMA30[i,res_term_K4_SMA30[i]]
}

adp_SMA30 <- ng.BHFDR(pval_SMA30)

res_term2_K4_SMA30 <- rep(NA,length=nrow(sma_K4))
for(i in 1:nrow(sma_K4)){
  if(adp_SMA30[i] < 0.001){
    res_term2_K4_SMA30[i] <- term[res_term_K4_SMA30[i]]
  }
}

table(res_term2_K4_SMA30)
plot(table(res_term2_K4_SMA30))
sum(table(res_term2_K4_SMA30))/nrow(sma_K4)


pdf("figs/climate_correlation/SMA7_Histone_cor_pval0.001_result.pdf",height=1.2,width=2)
set.panel(1,2)
par(mar=c(2.2,2,1,0.5))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)

names=c("Temperature","Precipitation","Day length","Sunlight hours")

barplot=barplot(table(res_term2_K27_SMA7)[c("temp","ppt","day","sun")],ylim=c(0,1500),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,1500,500), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-548,-548,-491,-615),labels=names,srt=45)
mtext(c("0","500","1,000","1,500"),at=seq(0,1500,500),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K27me3",side=3,line=0.1)

barplot=barplot(table(res_term2_K4_SMA7)[c("temp","ppt","day","sun")],ylim=c(0,1500),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,1500,500), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-548,-548,-491,-615),labels=names,srt=45)
mtext(c("0","500","1,000","1,500"),at=seq(0,1500,500),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K4me3",side=3,line=0.1)

dev.off()


pdf("figs/climate_correlation/SMA30_Histone_cor_pval0.001_result.pdf",height=1.2,width=2)
set.panel(1,2)
par(mar=c(2.2,2,1,0.5))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)

names=c("Temperature","Precipitation","Day length","Sunlight hours")

barplot=barplot(table(res_term2_K27_SMA30)[c("temp","ppt","day","sun")],ylim=c(0,1500),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,1500,500), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-548,-548,-491,-615),labels=names,srt=45)
mtext(c("0","500","1,000","1,500"),at=seq(0,1500,500),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K27me3",side=3,line=0.1)

barplot=barplot(table(res_term2_K4_SMA30)[c("temp","ppt","day","sun")],ylim=c(0,1500),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, label=F, at=barplot)
axis(side=2, at=seq(0,1500,500), pos=0, tcl=-0.1, labels=F)
text(c(-0.8,0.7,2.4,3.3),c(-548,-548,-491,-615),labels=names,srt=45)
mtext(c("0","500","1,000","1,500"),at=seq(0,1500,500),side=2,line=0.2,las=1)
mtext("No. of genes",side=2,line=1.3)
arrows(0, 0, 5.8, 0, length=0)
mtext("H3K4me3",side=3,line=0.1)

dev.off()
