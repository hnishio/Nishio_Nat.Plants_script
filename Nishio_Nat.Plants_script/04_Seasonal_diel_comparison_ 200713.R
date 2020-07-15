
library(LSD)
library(MASS)
library(fields)
library(exactRankTests)
library(lawstat)

setwd("../")
source("functions/ng.Colors.R")



##### Seasonal-diel comparison (all genes) #####
setwd("data")

K4.K27.seasonal <- read.csv("K4K27_data_stat.csv",header=T,sep=",")
K4.K27.diel <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")


# Amplitude of K4
dir.create("../figs/scatter/")

pdf("../figs/scatter/heatscatter_ampK4_K4K27_seasonal_diel_allgenes.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,5.5),ylim=c(0,5.5),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,5,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3","4","5")
mtext(side=1,at=seq(0,5,by=1),labelx,line=-0.35)
axis(2, at=seq(0,5,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,5,by=1),labelx,line=0.15,las=1)
mtext("H3K4me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.K27.diel$ampK4, K4.K27.seasonal$ampK4)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.6)
dev.off()


# Amplitude of K27
pdf("../figs/scatter/heatscatter_ampK27_K4K27_seasonal_diel_allgenes.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27.diel$ampK27, K4.K27.seasonal$ampK27, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,5.5),ylim=c(0,5.5), 
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,5,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3","4","5")
mtext(side=1,at=seq(0,5,by=1),labelx,line=-0.35)
axis(2, at=seq(0,5,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,5,by=1),labelx,line=0.15,las=1)
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
dir.create("../figs/boxplot")

pdf("../figs/boxplot/boxplot_ampK4_seasonal_diel_allgenes.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,2))
axis(1,at=1:2,labels=c("",""))
axis(2,las=1, labels=F)
b <-
boxplot(K4.K27.seasonal$ampK4, K4.K27.diel$ampK4, 
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08,2.5,-0.08,length=0,lwd=1)
arrows(0.42,-0.08,0.42,0.5,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.1,1,1.6,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.1,2,1.6,length=0,lwd=0.5)
arrows(1,1.6,2,1.6,length=0,lwd=0.5)
text(1.1,1.90,expression(paste(italic(P)," < ")),cex=0.8)
text(1.8,1.75,expression(paste("1.0 × ",10^-15)),cex=0.8)

text(c(0.17,1.7),c(-0.62,-0.4),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1","2"), side=2, line=0.15, at=c(0,1,2), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K4me3", side=3, line=0.1)
dev.off()


# amplitude of K27
pdf("../figs/boxplot/boxplot_ampK27_seasonal_diel_allgenes.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,las=1, labels=F)
b <-
boxplot(K4.K27.seasonal$ampK27, K4.K27.diel$ampK27, 
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08,2.5,-0.08,length=0,lwd=1)
arrows(0.42,-0.08,0.42,0.5,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.1,1,1.8,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.1,2,1.8,length=0,lwd=0.5)
arrows(1,1.8,2,1.8,length=0,lwd=0.5)
text(1.1,2.1,expression(paste(italic(P)," < ")),cex=0.8)
text(1.8,1.95,expression(paste("1.0 × ",10^-15)),cex=0.8)

text(c(0.17,1.7),c(-0.62,-0.4),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1","2"), side=2, line=0.15, at=c(0,1,2), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K27me3", side=3, line=0.1)
dev.off()



##### Amplitude ratio (seasonal / diel) #####
K4.S.D <- K4.K27.seasonal$ampK4 / K4.K27.diel$ampK4
K27.S.D <- K4.K27.seasonal$ampK27 / K4.K27.diel$ampK27

# Addition of the ratio to the last column
K4.S.enrich.diel.S.D <- cbind(K4.S.enrich.diel, K4.S.D)
K27.S.enrich.diel.S.D <- cbind(K27.S.enrich.diel, K27.S.D)

# Remove the genes with denominator = 0
K4.S.enrich.diel.S.D.rev <- 
  K4.S.enrich.diel.S.D[is.finite(K4.S.enrich.diel.S.D$K4.S.D),]
K27.S.enrich.diel.S.D.rev <- 
  K27.S.enrich.diel.S.D[is.finite(K27.S.enrich.diel.S.D$K27.S.D),]

# Ordering of genes according to the S-D ratio
sortlist <- order(K4.S.enrich.diel.S.D.rev$K4.S.D, decreasing=T)
K4.S.enrich.diel.S.D.rev.sort <- K4.S.enrich.diel.S.D.rev[sortlist,]
sortlist <- order(K27.S.enrich.diel.S.D.rev$K27.S.D, decreasing=T)
K27.S.enrich.diel.S.D.rev.sort <- K27.S.enrich.diel.S.D.rev[sortlist,]

write.csv(K4.S.enrich.diel.S.D.rev.sort, 
          file="K4K27_seasonal_diel_K4.S.enrich_SDratio.csv",quote=F, row.names=F)
write.csv(K27.S.enrich.diel.S.D.rev.sort, 
          file="K4K27_seasonal_diel_K27.S.enrich_SDratio.csv",quote=F, row.names=F)

K4.S.D <- K4.S.D[is.finite(K4.S.D)]
K27.S.D <- K27.S.D[is.finite(K27.S.D)]


# Brunner-Munzel test
brunner.munzel.test(x=K4.S.D, y=K27.S.D)


# Boxplot
pdf("../figs/boxplot/boxplot_SDratio_allgenes.pdf",width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,15))
axis(1,at=1:2,labels=c("",""))
axis(2,las=1, labels=F)
b <-
boxplot(K27.S.D,K4.S.D, las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.6,2.5,-0.6,length=0,lwd=1)
arrows(0.42,-0.6,0.42,1,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.6,1,15.5,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.6,2,15.5,length=0,lwd=0.5)
arrows(1,15.5,2,15.5,length=0,lwd=0.5)
text(1.3,16.7,expression(paste(italic(P)," < ","1.0 × ",10^-15)),cex=0.8)

text(c(0.1,1.2),c(-5,-4.6),c("H3K27me3", "H3K4me3"),srt=45)
mtext(c("0","5","10","15"), side=2, line=0.15, at=c(0,5,10,14.7), las=1)
mtext("Amplitude ratio", side=2, line=0.98)
mtext("(seasonal / diel)", side=2, line=0.55)
dev.off()



##### Seasonal-diel comparison (sesonally enriched genes) #####

# Merging seasonal and diel data for seasonally enriched genes
setwd("data")
K4.S.enrich <- read.csv("K4K27_data_stat_maxK4over2.csv",header=T,sep=",")
K27.S.enrich <- read.csv("K4K27_data_stat_maxK27over2.csv",header=T,sep=",")
K4.K27.diel <- read.csv("K4K27OMO48h_data_stat.csv",header=T,sep=",")

K4.S.enrich.diel <- 
  merge(K4.S.enrich, K4.K27.diel, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)

write.csv(K4.S.enrich.diel, 
          file="K4K27_seasonal_diel_K4.S.enrich.csv", 
          quote=F, row.names=F)

K27.S.enrich.diel <- 
  merge(K27.S.enrich, K4.K27.diel, 
        by=c("Ahal_ID","Araport11_ID","Reciprocal_flag","Locus_type","Symbol","Alias","Note"),
        all=F,sort=F)

write.csv(K27.S.enrich.diel, 
          file="K4K27_seasonal_diel_K27.S.enrich.csv", 
          quote=F, row.names=F)


# Heatscatter of seasonal and diel amplitudes
K4.S.enrich.diel <- 
  read.csv("K4K27_seasonal_diel_K4.S.enrich.csv",header=T,sep=",")
K27.S.enrich.diel <- 
  read.csv("K4K27_seasonal_diel_K27.S.enrich.csv",header=T,sep=",")

# Amplitude of K4
pdf("../figs/scatter/heatscatter_ampK4_K4K27_seasonal_diel_K4.S.enrich.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.S.enrich.diel$ampK4.y, K4.S.enrich.diel$ampK4.x, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,5.5),ylim=c(0,5.5),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,5,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3","4","5")
mtext(side=1,at=seq(0,5,by=1),labelx,line=-0.35)
axis(2, at=seq(0,5,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,5,by=1),labelx,line=0.15,las=1)
mtext("H3K4me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K4.S.enrich.diel$ampK4.y, K4.S.enrich.diel$ampK4.x)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.6)
dev.off()

# Amplitude of K27
pdf("../figs/scatter/heatscatter_ampK27_K4K27_seasonal_diel_K27.S.enrich.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K27.S.enrich.diel$ampK27.y, K27.S.enrich.diel$ampK27.x, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,5.5),ylim=c(0,5.5), 
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,5,by=1), label=F, tcl=-0.1)
labelx <- c("0","1","2","3","4","5")
mtext(side=1,at=seq(0,5,by=1),labelx,line=-0.35)
axis(2, at=seq(0,5,by=1), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,5,by=1),labelx,line=0.15,las=1)
mtext("H3K27me3",side=3,line=0)
mtext(expression(paste(log[2],"(diel amplitude)")),side=1,line=0.2)
mtext(expression(paste(log[2],"(seasonal amplitude)")),side=2,line=0.3)
d <- kde2d(K27.S.enrich.diel$ampK27.y, K27.S.enrich.diel$ampK27.x)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density",side=4,line=1.4)
dev.off()


# Wilcoxon rank sum test
wilcox.exact(x=K4.S.enrich.diel$ampK4.x, y=K4.S.enrich.diel$ampK4.y, 
             paired=F)
wilcox.exact(x=K27.S.enrich.diel$ampK27.x, y=K27.S.enrich.diel$ampK27.y, 
             paired=F)


# Brunner-Munzel test
brunner.munzel.test(x=K4.S.enrich.diel$ampK4.x, y=K4.S.enrich.diel$ampK4.y)
brunner.munzel.test(x=K27.S.enrich.diel$ampK27.x, y=K27.S.enrich.diel$ampK27.y)


# amplitude of K4
pdf("../figs/boxplot/boxplot_ampK4_K4K27_seasonal_diel_K4.S.enrich.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,2))
axis(1,at=1:2,labels=c("",""))
axis(2,las=1, labels=F)
b <-
boxplot(K4.S.enrich.diel$ampK4.x, K4.S.enrich.diel$ampK4.y, 
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08,2.5,-0.08,length=0,lwd=1)
arrows(0.42,-0.08,0.42,0.5,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.1,1,1.6,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.1,2,1.6,length=0,lwd=0.5)
arrows(1,1.6,2,1.6,length=0,lwd=0.5)
text(1.7,1.75,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62,-0.4),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1","2"), side=2, line=0.15, at=c(0,1,2), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K4me3", side=3, line=0)
dev.off()


# amplitude of K27
pdf("../figs/boxplot/boxplot_ampK27_K4K27_seasonal_diel_K27.S.enrich.pdf",
    width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,2))
names <- c("Seasonal", "Diel")
axis(1,at=1:2,labels=c("",""))
axis(2,las=1, labels=F)
b <-
boxplot(K27.S.enrich.diel$ampK27.x, K27.S.enrich.diel$ampK27.y, 
        las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.08,2.5,-0.08,length=0,lwd=1)
arrows(0.42,-0.08,0.42,0.5,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.1,1,1.8,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.1,2,1.8,length=0,lwd=0.5)
arrows(1,1.8,2,1.8,length=0,lwd=0.5)
text(1.7,1.95,expression(paste(italic(P)," < 0.001",)),cex=0.8)

text(c(0.17,1.7),c(-0.62,-0.4),c("Seasonal", "Diel"),srt=45)
mtext(c("0","1","2"), side=2, line=0.15, at=c(0,1,2), las=1)
mtext(expression(paste(log[2],"(amplitude)")), side=2, line=0.25)
mtext("H3K27me3", side=3, line=0)
dev.off()



##### Amplitude ratio (seasonal / diel) #####
K4.S.enrich.diel <- 
  read.csv("K4K27_seasonal_diel_K4.S.enrich.csv",header=T,sep=",")
K27.S.enrich.diel <- 
  read.csv("K4K27_seasonal_diel_K27.S.enrich.csv",header=T,sep=",")

# Ratio of seasonal / diel
K4.S.D <- K4.S.enrich.diel$ampK4.x / K4.S.enrich.diel$ampK4.y
K27.S.D <- K27.S.enrich.diel$ampK27.x / K27.S.enrich.diel$ampK27.y

# Addition of the ratio to the last column
K4.S.enrich.diel.S.D <- cbind(K4.S.enrich.diel, K4.S.D)
K27.S.enrich.diel.S.D <- cbind(K27.S.enrich.diel, K27.S.D)

# Remove the genes with denominator = 0
K4.S.enrich.diel.S.D.rev <- 
  K4.S.enrich.diel.S.D[is.finite(K4.S.enrich.diel.S.D$K4.S.D),]
K27.S.enrich.diel.S.D.rev <- 
  K27.S.enrich.diel.S.D[is.finite(K27.S.enrich.diel.S.D$K27.S.D),]

# Ordering of genes according to the S-D ratio
sortlist <- order(K4.S.enrich.diel.S.D.rev$K4.S.D, decreasing=T)
K4.S.enrich.diel.S.D.rev.sort <- K4.S.enrich.diel.S.D.rev[sortlist,]
sortlist <- order(K27.S.enrich.diel.S.D.rev$K27.S.D, decreasing=T)
K27.S.enrich.diel.S.D.rev.sort <- K27.S.enrich.diel.S.D.rev[sortlist,]

write.csv(K4.S.enrich.diel.S.D.rev.sort, 
          file="K4K27_seasonal_diel_K4.S.enrich_SDratio.csv",quote=F, row.names=F)
write.csv(K27.S.enrich.diel.S.D.rev.sort, 
          file="K4K27_seasonal_diel_K27.S.enrich_SDratio.csv",quote=F, row.names=F)


K4.S.D <- K4.S.D[is.finite(K4.S.D)]
K27.S.D <- K27.S.D[is.finite(K27.S.D)]


# Wilcoxon rank sum test
wilcox.exact(x=K4.S.D, y=K27.S.D, paired=F)


# Brunner-Munzel test
brunner.munzel.test(x=K4.S.D, y=K27.S.D)


# Boxplot
pdf("../figs/boxplot/boxplot_SDratio_S.enrich.pdf",width=0.9,height=1.3)
par(ps=6)
par(mgp=c(0,0.3,0))
par(mar=c(2.3,2.5,0.7,0.5))
par(tcl=-0.1)
par(xpd=T)
plot(0,0,type="n",xlab="",ylab="",axes=F,xlim=c(0.5,2.5),ylim=c(0,15))
axis(1,at=1:2,labels=c("",""))
axis(2,las=1, labels=F)
b <-
boxplot(K27.S.D,K4.S.D, las=1, xaxt="n",yaxt="n",add=T, cex=0.2, medlwd=2, outline=F, axes=F)
arrows(0.42,-0.6,2.5,-0.6,length=0,lwd=1)
arrows(0.42,-0.6,0.42,1,length=0,lwd=1)

arrows(1,b$stats[5,1]+0.6,1,15.5,length=0,lwd=0.5)
arrows(2,b$stats[5,2]+0.6,2,15.5,length=0,lwd=0.5)
arrows(1,15.5,2,15.5,length=0,lwd=0.5)
text(1.3,16.7,expression(paste(italic(P)," < ","1.0 × ",10^-15)),cex=0.8)

text(c(0.1,1.2),c(-5,-4.6),c("H3K27me3", "H3K4me3"),srt=45)
mtext(c("0","5","10","15"), side=2, line=0.15, at=c(0,5,10,14.7), las=1)
mtext("Amplitude ratio", side=2, line=0.98)
mtext("(seasonal / diel)", side=2, line=0.55)
dev.off()
