
library(LSD)
library(MASS)
library(fields)
library(exactRankTests)
library(lawstat)

source("functions/ng.Colors.R")



########## Sparse 1 ##########

##### Preparation for spline #####
setwd("data")

K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,c(11,14,17,20)],K4.K27.rep1[,c(11,14,17,20)],K4.K27.rep1[,c(11,14,17,20)])
K27.rep1 <- cbind(K4.K27.rep1[,c(26,29,32,35)],K4.K27.rep1[,c(26,29,32,35)],K4.K27.rep1[,c(26,29,32,35)])
K4.rep2 <- cbind(K4.K27.rep2[,c(11,14,17,20)],K4.K27.rep2[,c(11,14,17,20)],K4.K27.rep2[,c(11,14,17,20)])
K27.rep2 <- cbind(K4.K27.rep2[,c(26,29,32,35)],K4.K27.rep2[,c(26,29,32,35)],K4.K27.rep2[,c(26,29,32,35)])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- date[c(1,4,7,10)]
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
  spK4 <- smooth.spline(date,K4[i,1:24],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:24],spar=0.3)
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
write.csv(K4.K27, file="K4K27_data_stat_sparse1.csv", quote=F, row.names=F)



##### Seasonal diel comparison (all genes)
K4.K27.seasonal <- read.csv("K4K27_data_stat_sparse1.csv",header=T,sep=",")
K4.K27.diel <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")


# Amplitude of K4
pdf("../figs/scatter/heatscatter_ampK4_K4K27_seasonal_diel_allgenes_sparse1.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,3.5),ylim=c(0,3.5),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,3,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3")
mtext(side=1,at=seq(0,3,by=1),labelx,line=-0.35)
axis(2, at=seq(0,3,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,3,by=1),labelx,line=0.15,las=1)
mtext("H3K4me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.6)
dev.off()


# Amplitude of K27
pdf("../figs/scatter/heatscatter_ampK27_K4K27_seasonal_diel_allgenes_sparse1.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,3.5),ylim=c(0,3.5), 
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,3,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3")
mtext(side=1,at=seq(0,3,by=1),labelx,line=-0.35)
axis(2, at=seq(0,3,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,3,by=1),labelx,line=0.15,las=1)
mtext("H3K27me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.4)
dev.off()


# Brunner-Munzel test
brunner.munzel.test(x=K4.K27.seasonal$ampK4, y=K4.K27.diel$ampK4)
brunner.munzel.test(x=K4.K27.seasonal$ampK27, y=K4.K27.diel$ampK27)


# amplitude of K4
pdf("../figs/boxplot/boxplot_ampK4_seasonal_diel_allgenes_sparse1.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,1.2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,1,0.5),las=1,labels=F)
b <-
  boxplot(K4.K27.seasonal$ampK4, K4.K27.diel$ampK4, 
          las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08*1.2/2,2.5,-0.08*1.2/2,length=0,lwd=1)
arrows(0.42,-0.08*1.2/2,0.42,1.2,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.05,1,1.05,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.05,2,1.05,length=0,lwd=0.5)
arrows(1,1.05,2,1.05,length=0,lwd=0.5)
text(1.7,1.15,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62*1.2/2,-0.4*1.2/2),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1"), side=2, line=0.15, at=c(0,1), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K4me3", side=3, line=0)
dev.off()


# amplitude of K27
pdf("../figs/boxplot/boxplot_ampK27_seasonal_diel_allgenes_sparse1.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,1.2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,1,0.5),las=1,labels=F)
b <-
  boxplot(K4.K27.seasonal$ampK27, K4.K27.diel$ampK27, 
          las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08*1.2/2,2.5,-0.08*1.2/2,length=0,lwd=1)
arrows(0.42,-0.08*1.2/2,0.42,1.2,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.05,1,1.05,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.05,2,1.05,length=0,lwd=0.5)
arrows(1,1.05,2,1.05,length=0,lwd=0.5)
text(1.7,1.15,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62*1.2/2,-0.4*1.2/2),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1"), side=2, line=0.15, at=c(0,1), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K27me3", side=3, line=0)
dev.off()



##### Ratio of seasonal / diel #####
K4.S.D <- K4.K27.seasonal$ampK4 / K4.K27.diel$ampK4
K27.S.D <- K4.K27.seasonal$ampK27 / K4.K27.diel$ampK27

K4.S.D <- K4.S.D[is.finite(K4.S.D)]
K27.S.D <- K27.S.D[is.finite(K27.S.D)]


# Brunner-Munzel test
brunner.munzel.test(x=K4.S.D, y=K27.S.D)


# Boxplot
pdf("../figs/boxplot/boxplot_SDratio_allgenes_sparse1.pdf",width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,11))
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,10,5),las=1,labels=F)
b <-
  boxplot(K27.S.D,K4.S.D, las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.6*11/15,2.5,-0.6*11/15,length=0,lwd=1)
arrows(0.42,-0.6*11/15,0.42,11,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.4,1,11,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.4,2,11,length=0,lwd=0.5)
arrows(1,11,2,11,length=0,lwd=0.5)
text(1.3,12,expression(paste(italic(P)," < ","1.0 × ",10^-15)),cex=0.8)

text(c(0.1,1.2),c(-5*11/15,-4.6*11/15),c("H3K27me3", "H3K4me3"),srt=45)
mtext(c("0","5","10"), side=2, line=0.15, at=seq(0,10,5), las=1)
mtext("Amplitude ratio", side=2, line=0.98)
mtext("(seasonal / diel)", side=2, line=0.55)
dev.off()



########## Sparse 2 ##########

##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,c(12,15,18,21)],K4.K27.rep1[,c(12,15,18,21)],K4.K27.rep1[,c(12,15,18,21)])
K27.rep1 <- cbind(K4.K27.rep1[,c(27,30,33,36)],K4.K27.rep1[,c(27,30,33,36)],K4.K27.rep1[,c(27,30,33,36)])
K4.rep2 <- cbind(K4.K27.rep2[,c(12,15,18,21)],K4.K27.rep2[,c(12,15,18,21)],K4.K27.rep2[,c(12,15,18,21)])
K27.rep2 <- cbind(K4.K27.rep2[,c(27,30,33,36)],K4.K27.rep2[,c(27,30,33,36)],K4.K27.rep2[,c(27,30,33,36)])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- date[c(2,5,8,11)]
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
  spK4 <- smooth.spline(date,K4[i,1:24],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:24],spar=0.3)
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
write.csv(K4.K27, file="K4K27_data_stat_sparse2.csv", quote=F, row.names=F)



##### Seasonal diel comparison (all genes)
K4.K27.seasonal <- read.csv("K4K27_data_stat_sparse2.csv",header=T,sep=",")
K4.K27.diel <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")


# Amplitude of K4
pdf("../figs/scatter/heatscatter_ampK4_K4K27_seasonal_diel_allgenes_sparse2.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,3.5),ylim=c(0,3.5),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,3,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3")
mtext(side=1,at=seq(0,3,by=1),labelx,line=-0.35)
axis(2, at=seq(0,3,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,3,by=1),labelx,line=0.15,las=1)
mtext("H3K4me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.6)
dev.off()


# Amplitude of K27
pdf("../figs/scatter/heatscatter_ampK27_K4K27_seasonal_diel_allgenes_sparse2.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,3.5),ylim=c(0,3.5), 
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,3,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3")
mtext(side=1,at=seq(0,3,by=1),labelx,line=-0.35)
axis(2, at=seq(0,3,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,3,by=1),labelx,line=0.15,las=1)
mtext("H3K27me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.4)
dev.off()


# Brunner-Munzel test
brunner.munzel.test(x=K4.K27.seasonal$ampK4, y=K4.K27.diel$ampK4)
brunner.munzel.test(x=K4.K27.seasonal$ampK27, y=K4.K27.diel$ampK27)


# amplitude of K4
pdf("../figs/boxplot/boxplot_ampK4_seasonal_diel_allgenes_sparse2.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,1.2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,1,0.5),las=1,labels=F)
b <-
  boxplot(K4.K27.seasonal$ampK4, K4.K27.diel$ampK4, 
          las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08*1.2/2,2.5,-0.08*1.2/2,length=0,lwd=1)
arrows(0.42,-0.08*1.2/2,0.42,1.2,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.05,1,1.05,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.05,2,1.05,length=0,lwd=0.5)
arrows(1,1.05,2,1.05,length=0,lwd=0.5)
text(1.7,1.15,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62*1.2/2,-0.4*1.2/2),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1"), side=2, line=0.15, at=c(0,1), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K4me3", side=3, line=0)
dev.off()


# amplitude of K27
pdf("../figs/boxplot/boxplot_ampK27_seasonal_diel_allgenes_sparse2.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,1.2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,1,0.5),las=1,labels=F)
b <-
  boxplot(K4.K27.seasonal$ampK27, K4.K27.diel$ampK27, 
          las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08*1.2/2,2.5,-0.08*1.2/2,length=0,lwd=1)
arrows(0.42,-0.08*1.2/2,0.42,1.2,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.05,1,1.05,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.05,2,1.05,length=0,lwd=0.5)
arrows(1,1.05,2,1.05,length=0,lwd=0.5)
text(1.7,1.15,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62*1.2/2,-0.4*1.2/2),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1"), side=2, line=0.15, at=c(0,1), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K27me3", side=3, line=0)
dev.off()



##### Ratio of seasonal / diel #####
K4.S.D <- K4.K27.seasonal$ampK4 / K4.K27.diel$ampK4
K27.S.D <- K4.K27.seasonal$ampK27 / K4.K27.diel$ampK27

K4.S.D <- K4.S.D[is.finite(K4.S.D)]
K27.S.D <- K27.S.D[is.finite(K27.S.D)]


# Brunner-Munzel test
brunner.munzel.test(x=K4.S.D, y=K27.S.D)


# Boxplot
pdf("../figs/boxplot/boxplot_SDratio_allgenes_sparse2.pdf",width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,11))
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,10,5),las=1,labels=F)
b <-
  boxplot(K27.S.D,K4.S.D, las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.6*11/15,2.5,-0.6*11/15,length=0,lwd=1)
arrows(0.42,-0.6*11/15,0.42,11,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.4,1,11,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.4,2,11,length=0,lwd=0.5)
arrows(1,11,2,11,length=0,lwd=0.5)
text(1.3,12,expression(paste(italic(P)," < ","1.0 × ",10^-15)),cex=0.8)

text(c(0.1,1.2),c(-5*11/15,-4.6*11/15),c("H3K27me3", "H3K4me3"),srt=45)
mtext(c("0","5","10"), side=2, line=0.15, at=seq(0,10,5), las=1)
mtext("Amplitude ratio", side=2, line=0.98)
mtext("(seasonal / diel)", side=2, line=0.55)
dev.off()



########## Sparse 3 ##########

##### Preparation for spline #####
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,c(13,16,19,22)],K4.K27.rep1[,c(13,16,19,22)],K4.K27.rep1[,c(13,16,19,22)])
K27.rep1 <- cbind(K4.K27.rep1[,c(28,31,34,37)],K4.K27.rep1[,c(28,31,34,37)],K4.K27.rep1[,c(28,31,34,37)])
K4.rep2 <- cbind(K4.K27.rep2[,c(13,16,19,22)],K4.K27.rep2[,c(13,16,19,22)],K4.K27.rep2[,c(13,16,19,22)])
K27.rep2 <- cbind(K4.K27.rep2[,c(28,31,34,37)],K4.K27.rep2[,c(28,31,34,37)],K4.K27.rep2[,c(28,31,34,37)])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- date[c(3,6,9,12)]
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
  spK4 <- smooth.spline(date,K4[i,1:24],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:24],spar=0.3)
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
write.csv(K4.K27, file="K4K27_data_stat_sparse3.csv", quote=F, row.names=F)



##### Seasonal diel comparison (all genes)
K4.K27.seasonal <- read.csv("K4K27_data_stat_sparse3.csv",header=T,sep=",")
K4.K27.diel <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")


# Amplitude of K4
pdf("../figs/scatter/heatscatter_ampK4_K4K27_seasonal_diel_allgenes_sparse3.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,3.5),ylim=c(0,3.5),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,3,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3")
mtext(side=1,at=seq(0,3,by=1),labelx,line=-0.35)
axis(2, at=seq(0,3,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,3,by=1),labelx,line=0.15,las=1)
mtext("H3K4me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.6)
dev.off()


# Amplitude of K27
pdf("../figs/scatter/heatscatter_ampK27_K4K27_seasonal_diel_allgenes_sparse3.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,3.5),ylim=c(0,3.5), 
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,3,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3")
mtext(side=1,at=seq(0,3,by=1),labelx,line=-0.35)
axis(2, at=seq(0,3,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,3,by=1),labelx,line=0.15,las=1)
mtext("H3K27me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.4)
dev.off()


# Brunner-Munzel test
brunner.munzel.test(x=K4.K27.seasonal$ampK4, y=K4.K27.diel$ampK4)
brunner.munzel.test(x=K4.K27.seasonal$ampK27, y=K4.K27.diel$ampK27)


# amplitude of K4
pdf("../figs/boxplot/boxplot_ampK4_seasonal_diel_allgenes_sparse3.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,1.2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,1,0.5),las=1,labels=F)
b <-
  boxplot(K4.K27.seasonal$ampK4, K4.K27.diel$ampK4, 
          las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08*1.2/2,2.5,-0.08*1.2/2,length=0,lwd=1)
arrows(0.42,-0.08*1.2/2,0.42,1.2,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.05,1,1.05,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.05,2,1.05,length=0,lwd=0.5)
arrows(1,1.05,2,1.05,length=0,lwd=0.5)
text(1.7,1.15,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62*1.2/2,-0.4*1.2/2),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1"), side=2, line=0.15, at=c(0,1), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K4me3", side=3, line=0)
dev.off()


# amplitude of K27
pdf("../figs/boxplot/boxplot_ampK27_seasonal_diel_allgenes_sparse3.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,1.2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,1,0.5),las=1,labels=F)
b <-
  boxplot(K4.K27.seasonal$ampK27, K4.K27.diel$ampK27, 
          las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08*1.2/2,2.5,-0.08*1.2/2,length=0,lwd=1)
arrows(0.42,-0.08*1.2/2,0.42,1.2,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.05,1,1.05,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.05,2,1.05,length=0,lwd=0.5)
arrows(1,1.05,2,1.05,length=0,lwd=0.5)
text(1.7,1.15,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62*1.2/2,-0.4*1.2/2),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1"), side=2, line=0.15, at=c(0,1), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K27me3", side=3, line=0)
dev.off()



##### Ratio of seasonal / diel #####
K4.S.D <- K4.K27.seasonal$ampK4 / K4.K27.diel$ampK4
K27.S.D <- K4.K27.seasonal$ampK27 / K4.K27.diel$ampK27

K4.S.D <- K4.S.D[is.finite(K4.S.D)]
K27.S.D <- K27.S.D[is.finite(K27.S.D)]


# Brunner-Munzel test
brunner.munzel.test(x=K4.S.D, y=K27.S.D)


# Boxplot
pdf("../figs/boxplot/boxplot_SDratio_allgenes_sparse3.pdf",width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,11))
axis(1,at=1:2,labels=c("",""))
axis(2,at=seq(0,10,5),las=1,labels=F)
b <-
  boxplot(K27.S.D,K4.S.D, las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.6*11/15,2.5,-0.6*11/15,length=0,lwd=1)
arrows(0.42,-0.6*11/15,0.42,11,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.4,1,11,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.4,2,11,length=0,lwd=0.5)
arrows(1,11,2,11,length=0,lwd=0.5)
text(1.3,12,expression(paste(italic(P)," < ","1.0 × ",10^-15)),cex=0.8)

text(c(0.1,1.2),c(-5*11/15,-4.6*11/15),c("H3K27me3", "H3K4me3"),srt=45)
mtext(c("0","5","10"), side=2, line=0.15, at=seq(0,10,5), las=1)
mtext("Amplitude ratio", side=2, line=0.98)
mtext("(seasonal / diel)", side=2, line=0.55)
dev.off()
