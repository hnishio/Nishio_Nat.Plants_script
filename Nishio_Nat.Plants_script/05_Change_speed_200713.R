
library(LSD)
library(MASS)
library(fields)
library(exactRankTests)

setwd("../")
source("functions/ng.Colors.R")



##### Illustration of increase and decrease durations #####
setwd("data")

x<-seq(-pi/2,pi*3/2,pi/334)
s<-sin(x)

pdf("../figs/others/Inc_dec_duration_illustration.pdf",height=1.1,width=1.7)
par(mar=c(1,2,0.9,1.1))
par(mgp=c(3,0.15,0.15))
par(ps=6)
par(cex=1)

plot(s,type="l",xlim=c(0,668),ylim=c(-1.2,1.2),
     xlab="",ylab="",ann=F,axes=F)
day<-c(365,31,28,31,30,31,30,31,31,30,31,30,31)
par(new=T)
plot(NULL, xlim=c(365,730),ylim=c(0,1.2),xlab="",ylab="",ann=F,axes=F)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
       mean(place[9:10]),mean(place[11:12]))
labelx<-c("1","3","5","7","9","11")
axis(side=1,at=place,label=F,tcl=-0.1,pos=0)
mtext(side=1,at=lit,labelx,line=-0.5)
axis(side=2,at=c(0.1,0.2,1,1.1),tcl=-0.1,las=2,pos=365,labels=F)
mtext(c("0%","10%","90%","100%"),at=c(0.08,0.22,0.98,1.12),side=2,line=-0.05,las=1)
mtext("Month",side=1,line=-0.1)
mtext("Modification level",side=2,line=0.9)
arrows(365,1.2,730,1.2,length=0)
arrows(730,-0.2,730,1.2,length=0)
arrows(365,-0.2,365,1.2,length=0)
arrows(365,1,730,1,length=0,col="black",lwd=0.5,lty="dashed")
arrows(365,0.2,730,0.2,length=0,col="black",lwd=0.5,lty="dashed")

	p <- (asin(0.8)+pi/2)*365/(pi*2)
	p2 <- 365+p
	p3 <- 730-p
	p <- (asin(-0.8)+pi/2)*365/(pi*2)
	p1 <- 365+p
	p4 <- 730-p	
	arrows(p1,0,p1,1.2,length=0,col="black",lwd=0.5,lty="dashed")
	arrows(p2,0,p2,1.2,length=0,col="black",lwd=0.5,lty="dashed")
	arrows(p3,0,p3,1.2,length=0,col="black",lwd=0.5,lty="dashed")
	arrows(p4,0,p4,1.2,length=0,col="black",lwd=0.5,lty="dashed")

par(xpd=T)
arrows(p1,0.85,p2,0.85,length=0.02,lwd=1,code=3,col="orange")
arrows((p1+p2)/2,0.85,(p1+p2)/2,1.3,length=0,lwd=0.3,col="orange")
arrows(p3,0.85,p4,0.85,length=0.02,lwd=1,code=3,col="orange")
arrows((p3+p4)/2,0.85,(p3+p4)/2,1.3,length=0,lwd=0.3,col="orange")
text((p1+p2)/2,1.5,"Increase",col="orange")
text((p1+p2)/2,1.37,"duration",col="orange")
text((p3+p4)/2,1.5,"Decrease",col="orange")
text((p3+p4)/2,1.37,"duration",col="orange")
dev.off()



###### Integration of H3K4me3 and H3K27me3 data #####

# rep1
K4 <- read.csv("K4_rep1_rpkm_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_rep1_rpkm_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_rep1_rpkm_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)

#rep2
K4 <- read.csv("K4_rep2_rpkm_annotation.csv", sep=",", header=T)
K27 <- read.csv("K27_rep2_rpkm_annotation.csv", sep=",", header=T)
K4.K27 <- merge(K4, K27, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                              "Locus_type","Symbol","Alias","Note"),all=F,sort=F)
write.csv(K4.K27, 
          file="K4K27_rep2_rpkm_annotation.csv", 
          quote=F, row.names=F)
nrow(K4.K27)



##### Increase and decrease duration of K4 and K27 #####
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm_annotation.csv",
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

risetimeK4 <- rep(NA,length=nrow(K4))
falltimeK4 <- rep(NA,length=nrow(K4))
risetimeK27 <- rep(NA,length=nrow(K4))
falltimeK27 <- rep(NA,length=nrow(K4))

K4.K27 <- cbind(title, 
                risetimeK4, falltimeK4,
                risetimeK27, falltimeK27)


1:nrow(K4.K27)
for(i in 19146){
	
  spK4 <- smooth.spline(date,K4[i,1:72],spar=0.3)
  spK27 <- smooth.spline(date,K27[i,1:72],spar=0.3)
  length <- 1000
  x <- seq(365,730,length=length)
  predK4 <- predict(spK4,x)
  predK27 <- predict(spK27,x)


# rise and fall of K4  
  rangeK4 <- max(predK4$y)-min(predK4$y)  
  maxK4time <- 
    predK4$x[abs(predK4$y - (max(predK4$y)-rangeK4*0.1)) < rangeK4*0.01]
  minK4time <-   
    predK4$x[abs(predK4$y - (min(predK4$y)+rangeK4*0.1)) < rangeK4*0.01]
 
  if(max(minK4time) < min(maxK4time)){
  	risetimeK4 <- min(maxK4time) - max(minK4time)
  	falltimeK4 <- (730-max(maxK4time))+(min(minK4time)-365)
  }else if(max(maxK4time) < min(minK4time)){
  	risetimeK4 <- (730-max(minK4time))+(min(maxK4time)-365)
  	falltimeK4 <- min(minK4time) - max(maxK4time)  	
  }else if((min(maxK4time) < min(minK4time)) && (max(minK4time) < max(maxK4time))){
    risetimeK4 <- min(maxK4time[maxK4time > max(minK4time)]) - max(minK4time)
    falltimeK4 <- min(minK4time) - max(maxK4time[maxK4time < min(minK4time)]) 
  }else if((min(minK4time) < min(maxK4time)) && (max(maxK4time) < max(minK4time))){
    risetimeK4 <- min(maxK4time) - max(minK4time[minK4time < min(maxK4time)])
    falltimeK4 <- min(minK4time[minK4time > max(maxK4time)]) - max(maxK4time)
  }else{
  	risetimeK4 <- NA
  	falltimeK4 <- NA
  }

K4.K27$risetimeK4[i] <- risetimeK4
K4.K27$falltimeK4[i] <- falltimeK4
  

# rise and fall of K27  
  rangeK27 <- max(predK27$y)-min(predK27$y)  
  maxK27time <- 
    predK27$x[abs(predK27$y - (max(predK27$y)-rangeK27*0.1)) < rangeK27*0.01]
  minK27time <-   
    predK27$x[abs(predK27$y - (min(predK27$y)+rangeK27*0.1)) < rangeK27*0.01]
 
  if(max(minK27time) < min(maxK27time)){
  	risetimeK27 <- min(maxK27time) - max(minK27time)
  	falltimeK27 <- (730-max(maxK27time))+(min(minK27time)-365)
  }else if(max(maxK27time) < min(minK27time)){
  	risetimeK27 <- (730-max(minK27time))+(min(maxK27time)-365)
  	falltimeK27 <- min(minK27time) - max(maxK27time)  	
  }else if((min(maxK27time) < min(minK27time)) && (max(minK27time) < max(maxK27time)) &&
  (length(maxK27time[min(minK27time) < maxK27time && maxK27time < max(minK27time)])==0)){
    risetimeK27 <- min(maxK27time[maxK27time > max(minK27time)]) - max(minK27time)
    falltimeK27 <- min(minK27time) - max(maxK27time[maxK27time < min(minK27time)]) 
  }else if((min(minK27time) < min(maxK27time)) && (max(maxK27time) < max(minK27time)) &&
  (length(minK27time[min(maxK27time) < minK27time && minK27time < max(maxK27time)])==0)){
    risetimeK27 <- min(maxK27time) - max(minK27time[minK27time < min(maxK27time)])
    falltimeK27 <- min(minK27time[minK27time > max(maxK27time)]) - max(maxK27time)
  }else{
  	risetimeK27 <- NA
  	falltimeK27 <- NA
  }

K4.K27$risetimeK27[i] <- risetimeK27
K4.K27$falltimeK27[i] <- falltimeK27

}

write.csv(K4.K27, file="risefalltimeK4K27.csv", quote=F, row.names=F)


# Confirmation of graph
plot(predK4)
arrows(0,(max(predK4$y)-rangeK4*0.1),800,(max(predK4$y)-rangeK4*0.1),col="red")
arrows(0,(min(predK4$y)+rangeK4*0.1),800,(min(predK4$y)+rangeK4*0.1),col="red")

plot(predK27)
arrows(0,(max(predK27$y)-rangeK27*0.1),800,(max(predK27$y)-rangeK27*0.1),col="red")
arrows(0,(min(predK27$y)+rangeK27*0.1),800,(min(predK27$y)+rangeK27*0.1),col="red")



##### Heatscatter of the time #####
K4.K27 <- read.csv("risefalltimeK4K27.csv",header=T,sep=",")
K4.K27[which(K4.K27==Inf, TRUE)] <- NA
K4.K27 <- na.omit(K4.K27[,8:11])


# K4
pdf("../figs/scatter/heatscatter_inc_dec_duration_K4.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27$risetimeK4, K4.K27$falltimeK4, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,315),ylim=c(0,315),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,300,by=100), label=F, tcl=-0.1)
labelx <- c("0","100","200","300")
mtext(side=1,at=seq(0,300,by=100),labelx,line=-0.35)
axis(2, at=seq(0,300,by=100), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,300,by=100),labelx,line=0.15,las=1)
mtext("Increase duration (days)",side=1,line=0.15)
mtext("Decrease duration (days)",side=2,line=0.8)
mtext("H3K4me3",side=3,line=0)
d <- kde2d(K4.K27$risetimeK4, K4.K27$falltimeK4)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density of genes",side=4,line=2.1)
dev.off()


# K27
pdf("../figs/scatter/heatscatter_inc_dec_duration_K27.pdf",
    width=1.75,height=1.7)
par(oma=c(0,0,0,0))
par(mar=c(2.4,1.8,2,3))
par(ps=6)
par(mgp=c(0,0.3,0))
par(xpd=T)
heatscatter(K4.K27$risetimeK27, K4.K27$falltimeK27, 
            pch=".", las=1, tcl=-0.1, xlab="",ylab="",xlim=c(0,315),ylim=c(0,315),
            main="",xaxt="n",yaxt="n",colpal=ng.po.colors(64))
axis(1, at=seq(0,300,by=100), label=F, tcl=-0.1)
labelx <- c("0","100","200","300")
mtext(side=1,at=seq(0,300,by=100),labelx,line=-0.35)
axis(2, at=seq(0,300,by=100), label=F, tcl=-0.1)
mtext(side=2,at=seq(0,300,by=100),labelx,line=0.15,las=1)
mtext("Increase duration (days)",side=1,line=0.15)
mtext("Decrease duration (days)",side=2,line=0.8)
mtext("H3K27me3",side=3,line=0)
d <- kde2d(K4.K27$risetimeK27, K4.K27$falltimeK27)
image.plot(legend.only=T,zlim=c(min(d$z),max(d$z)),legend.width=0.3,tcl=-0.1,
           col=ng.po.colors(64),legend.mar=3.4)
mtext("Density of genes",side=4,line=2.1)
dev.off()
