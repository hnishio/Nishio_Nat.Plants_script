
library(season)
library(VennDiagram)
library(fields)



###### Integration of H3K4me3 and H3K27me3 data #####
setwd("data")

# Preparation for cosinor
K4.K27.rep1 <- 
  read.csv("K4K27_rep1_rpkm+1_log2_annotation.csv",
           header=T,sep=",")
K4.K27.rep2 <- 
  read.csv("K4K27_rep2_rpkm+1_log2_annotation.csv",
           header=T,sep=",")

K4.K27.duplicate <- merge(K4.K27.rep1, K4.K27.rep2, by=c("Ahal_ID","Araport11_ID","Reciprocal_flag",
                                                         "Locus_type","Symbol","Alias","Note"),all=F,sort=F)

K4.rep1 <- cbind(K4.K27.rep1[,11:22])
K27.rep1 <- cbind(K4.K27.rep1[,26:37])
K4.rep2 <- cbind(K4.K27.rep2[,11:22])
K27.rep2 <- cbind(K4.K27.rep2[,26:37])

K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c("2012-11-06", "2012-12-04", "2013-01-08", "2013-02-05", "2013-03-05", "2013-04-02",
          "2013-04-30", "2013-05-28", "2013-07-02", "2013-07-30", "2013-08-27", "2013-09-24")
date <- as.Date(c(date,date))
title <- K4.K27.duplicate[,1:7]


K4.mat <- matrix(nrow=nrow(K4),ncol=9)
K27.mat <- matrix(nrow=nrow(K4),ncol=9)


for(i in 1:nrow(K4)){
  his.data <- data.frame(date=date,K4=K4[i,],K27=K27[i,])
  
  res.K4 <- cosinor(K4~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K4)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
 
  lphase <- gsub("Month = ", "2012-", summary(res.K4)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K4.mat[i,] <- c(summary(res.K4)$amp*2,peak,trough,
                  summary(res.K4$glm)$coefficients[,"Estimate"],
                  summary(res.K4$glm)$coefficients[,"Pr(>|t|)"])
  
  res.K27 <- cosinor(K27~1,date="date",data=his.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.K27)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.K27)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  K27.mat[i,] <- c(summary(res.K27)$amp*2,peak,trough,
                  summary(res.K27$glm)$coefficients[,"Estimate"],
                  summary(res.K27$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(K4.mat) <- c("amp.K4","peak.K4","trough.K4",
                      "coef.inter.K4","coef.cosw.K4","coef.sinw.K4",
                      "pval.inter.K4","pval.cosw.K4","pval.sinw.K4")
colnames(K27.mat) <- c("amp.K27","peak.K27","trough.K27",
                      "coef.inter.K27","coef.cosw.K27","coef.sinw.K27",
                      "pval.inter.K27","pval.cosw.K27","pval.sinw.K27")

K4.K27.mat <- cbind(title,K4.mat,K27.mat)
write.csv(K4.K27.mat, file="K4K27_cosinor.csv", quote=F, row.names=F)


# Venn diagram
dir.create("../figs/cosinor_venn")

K4.K27.mat <- read.csv("K4K27_cosinor.csv",header=T,sep=",")
K4.CPMover1 <- read.csv("seasonalK4_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.CPMover1 <- read.csv("seasonalK27_CPMover1_edgeR_result.csv",header=T,sep=",")
K4 <- merge(K4.K27.mat,K4.CPMover1,by="Ahal_ID",all=F,sort=F)
K27 <- merge(K4.K27.mat,K27.CPMover1,by="Ahal_ID",all=F,sort=F)
K4.K27 <- merge(K4,K27,by="Ahal_ID",all=F,sort=F)
write.csv(K4[,1:ncol(K4.K27.mat)], file="K4K27_seasonal_cosinor_K4_CPMover1.csv", quote=F, row.names=F)
write.csv(K27[,1:ncol(K4.K27.mat)], file="K4K27_seasonal_cosinor_K27_CPMover1.csv", quote=F, row.names=F)

K4_season_osci <- length(which(K4$pval.cosw.K4<0.025 | K4$pval.sinw.K4<0.025))
K27_season_osci <- length(which(K27$pval.cosw.K27<0.025 | K27$pval.sinw.K27<0.025))
K4K27_season_osci <- length(which((K4.K27$pval.cosw.K4.x<0.025 | K4.K27$pval.sinw.K4.x<0.025) &
                                    (K4.K27$pval.cosw.K27.x<0.025 | K4.K27$pval.sinw.K27.x<0.025)))

pdf("../figs/cosinor_venn/venn_K4K27_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=K4_season_osci, area2=K27_season_osci, cross.area=K4K27_season_osci, 
  cex=1, category=c("K4 seasonal","K27 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()



##### Diel oscillation #####
# Load the diel data
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

K4.K27.rep1 <- K4.K27.rep1[order(K4.K27.rep1$Ahal_ID),]
K4.K27.rep2 <- K4.K27.rep2[order(K4.K27.rep2$Ahal_ID),]
K4.K27.rep3 <- K4.K27.rep3[order(K4.K27.rep3$Ahal_ID),]
K4.K27.rep4 <- K4.K27.rep4[order(K4.K27.rep4$Ahal_ID),]

K4.K27.quadruplicate <- cbind(K4.K27.rep1,K4.K27.rep2[,8:29],
                              K4.K27.rep3[,8:29],K4.K27.rep4[,8:29])

K4.rep1 <- cbind(K4.K27.rep1[,11:18])
K27.rep1 <- cbind(K4.K27.rep1[,22:29])
K4.rep2 <- cbind(K4.K27.rep2[,11:18])
K27.rep2 <- cbind(K4.K27.rep2[,22:29])
K4.rep3 <- cbind(K4.K27.rep3[,11:18])
K27.rep3 <- cbind(K4.K27.rep3[,22:29])
K4.rep4 <- cbind(K4.K27.rep4[,11:18])
K27.rep4 <- cbind(K4.K27.rep4[,22:29])

K4.diel <- cbind(K4.rep1,K4.rep2,K4.rep3,K4.rep4)
K27.diel <- cbind(K27.rep1,K27.rep2,K27.rep3,K27.rep4)
K4.diel <- as.matrix(K4.diel)
K27.diel <- as.matrix(K27.diel)
time <- c("2016-09-13 18:00", "2016-09-13 24:00", "2016-09-14 6:00", "2016-09-14 12:00", "2016-09-14 18:00", "2016-09-14 24:00",
          "2016-09-15 6:00", "2016-09-15 12:00")
time <- as.POSIXct(c(time,time,time,time))
title <- K4.K27.quadruplicate[,1:7]


K4.diel.mat <- matrix(nrow=nrow(K4.diel),ncol=9)
K27.diel.mat <- matrix(nrow=nrow(K4.diel),ncol=9)


for(i in 1:nrow(K4.diel)){
  his.data <- data.frame(time=time,K4.diel=K4.diel[i,],K27.diel=K27.diel[i,])
  
  res.K4.diel <- cosinor(K4.diel~1,date="time",type="hourly",data=his.data)
  
  peak <- as.numeric(gsub("Hour = ", "", summary(res.K4.diel)$phase))
  trough <- as.numeric(gsub("Hour = ", "", summary(res.K4.diel)$lphase))
  
  K4.diel.mat[i,] <- c(summary(res.K4.diel)$amp*2,peak,trough,
                  summary(res.K4.diel$glm)$coefficients[,"Estimate"],
                  summary(res.K4.diel$glm)$coefficients[,"Pr(>|t|)"])
  
  res.K27.diel <- cosinor(K27.diel~1,date="time",type="hourly",data=his.data)
  
  peak <- as.numeric(gsub("Hour = ", "", summary(res.K27.diel)$phase))
  trough <- as.numeric(gsub("Hour = ", "", summary(res.K27.diel)$lphase))
  
  K27.diel.mat[i,] <- c(summary(res.K27.diel)$amp*2,peak,trough,
                       summary(res.K27.diel$glm)$coefficients[,"Estimate"],
                       summary(res.K27.diel$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(K4.diel.mat) <- c("amp.K4.diel","peak.K4.diel","trough.K4.diel",
                      "coef.inter.K4.diel","coef.cosw.K4.diel","coef.sinw.K4.diel",
                      "pval.inter.K4.diel","pval.cosw.K4.diel","pval.sinw.K4.diel")
colnames(K27.diel.mat) <- c("amp.K27.diel","peak.K27.diel","trough.K27.diel",
                       "coef.inter.K27.diel","coef.cosw.K27.diel","coef.sinw.K27.diel",
                       "pval.inter.K27.diel","pval.cosw.K27.diel","pval.sinw.K27.diel")

K4.diel.K27.diel.mat <- cbind(title,K4.diel.mat,K27.diel.mat)
write.csv(K4.diel.K27.diel.mat, file="K4K27_diel_cosinor.csv", quote=F, row.names=F)


# Venn diagram
K4.K27.mat <- read.csv("K4K27_diel_cosinor.csv",header=T,sep=",")
K4.CPMover1 <- read.csv("dielK4_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.CPMover1 <- read.csv("dielK27_CPMover1_edgeR_result.csv",header=T,sep=",")
K4 <- merge(K4.K27.mat,K4.CPMover1,by="Ahal_ID",all=F,sort=F)
K27 <- merge(K4.K27.mat,K27.CPMover1,by="Ahal_ID",all=F,sort=F)
K4.K27 <- merge(K4,K27,by="Ahal_ID",all=F,sort=F)
write.csv(K4[,1:ncol(K4.K27.mat)], file="K4K27_diel_cosinor_K4_CPMover1.csv", quote=F, row.names=F)
write.csv(K27[,1:ncol(K4.K27.mat)], file="K4K27_diel_cosinor_K27_CPMover1.csv", quote=F, row.names=F)

K4_diel_osci <- length(which(K4$pval.cosw.K4.diel<0.025 | K4$pval.sinw.K4.diel<0.025))
K27_diel_osci <- length(which(K27$pval.cosw.K27.diel<0.025 | K27$pval.sinw.K27.diel<0.025))
K4K27_diel_osci <- length(which((K4.K27$pval.cosw.K4.diel.x<0.025 | K4.K27$pval.sinw.K4.diel.x<0.025) &
                                    (K4.K27$pval.cosw.K27.diel.x<0.025 | K4.K27$pval.sinw.K27.diel.x<0.025)))


pdf("../figs/cosinor_venn/venn_K4K27_diel_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=K4_diel_osci, area2=K27_diel_osci, cross.area=K4K27_diel_osci, 
  cex=1, category=c("K4 diel","K27 diel"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()



##### Graphs of seasonal and diel change #####
dir.create("../figs/cosinor_seasonal_diel_change")

day<-c(0,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]))


for(i in c(1263,1268,1821,2309,7214,8236,9640,9784,9826,10555,10811,13729,14931,15060,16823,17199,17215,17335,18162,18589,19190,20174,20212,20684,27591,28191)){
    pdf(paste("../figs/cosinor_seasonal_diel_change/cosinorK27K4seasonal_diel_", 
            title$Ahal_ID[i], "_", title$Symbol[i], ".pdf", sep=""),
      height=0.45,width=3.1) 
  set.panel(1,4,relax=T)
  par(oma=c(0,0.1,0,0))
  par(mar=c(0.6,0.5,0.2,0.1))
  par(mgp=c(3,0.13,0.15))
  par(ps=6)
  par(cex=1)
  
  # Cosinor model for seasonal change
  his.data <- data.frame(date=date,K4=K4[i,],K27=K27[i,])
  res.K4 <- cosinor(K4~1,date="date",data=his.data)
  res.K27 <- cosinor(K27~1,date="date",data=his.data)
  
  date2 <- (date - as.Date("2012-10-31"))[1:12]
  date3 <- seq(1,365,by=1)
  ft <- (date3-1)/365
  omega.t <- 2*pi*ft
  
  y.K4 <- summary(res.K4$glm)$coefficients[2,"Estimate"]*cos(omega.t) + 
    summary(res.K4$glm)$coefficients[3,"Estimate"]*sin(omega.t)
  y.K27 <- summary(res.K27$glm)$coefficients[2,"Estimate"]*cos(omega.t) + 
    summary(res.K27$glm)$coefficients[3,"Estimate"]*sin(omega.t)
  
  
  # Cosinor model for diel change
  his.data <- data.frame(time=time,K4.diel=K4.diel[i,],K27.diel=K27.diel[i,])
  res.K4.diel <- cosinor(K4.diel~1,date="time",type="hourly",data=his.data)
  res.K27.diel <- cosinor(K27.diel~1,date="time",type="hourly",data=his.data)
  
  time2 <- (time-as.POSIXct("2016-09-13 00:00"))
  time3 <- seq(as.numeric(time2[1]),as.numeric(time2[8]),by=1)
  ft <- (time3-1)/24
  omega.t <- 2*pi*ft
  
  y.K4.diel <- summary(res.K4.diel$glm)$coefficients[2,"Estimate"]*cos(omega.t) + 
    summary(res.K4.diel$glm)$coefficients[3,"Estimate"]*sin(omega.t)
  y.K27.diel <- summary(res.K27.diel$glm)$coefficients[2,"Estimate"]*cos(omega.t) + 
    summary(res.K27.diel$glm)$coefficients[3,"Estimate"]*sin(omega.t)
  
  
  # K27 plot (seasonal)
  lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
         mean(place[9:10]),mean(place[11:12]))
  labelx<-c("11","1","3","5","7","9")
  
  top <- 3
  bottom <- -3
  
  plot(rep(date2,2),scale(K27[i,]),
       type="p",las=1,tcl=-0.2,xlim=c(0,334),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  lines(c(1:334), scale(y.K27[c(305:365,1:273)]), col=rgb(0,0,0.7))
  
  axis(side=1,at=place,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  axis(side=2,at=c(-3,0,3),tcl=-0.1,las=2,pos=0)
  arrows(0,bottom,0,top,length=0)	
  arrows(0,top,334,top,length=0)
  arrows(334,bottom,334,top,length=0)	
  
  
  # K27 plot (diel)
  tp<-c(12,24,48,63)
  lit<-c(mean(tp[1:2]),mean(tp[2:3]),mean(tp[3:4]))
  labelx<-c("13","14","15")

  plot(scale(K27.diel[i,]),type="n",las=1,tcl=-0.2,xlim=c(12,63),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.5)
  rect(18.17,bottom,29.7,top,col="gray85",border="transparent")
  rect(42.15,bottom,53.72,top,col="gray85",border="transparent")
  par(new=T)
  plot(time2,scale(K27.diel[i,]),
       type="p",las=1,tcl=-0.2,xlim=c(12,63),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  lines(time3, scale(y.K27.diel), col=rgb(0,0,0.7))
  
  axis(side=1,at=tp,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  axis(side=2,at=c(-3,0,3),tcl=-0.1,las=2,pos=12)
  arrows(12,bottom,12,top,length=0)	
  arrows(12,top,63,top,length=0)
  arrows(63,bottom,63,top,length=0)	
  
  
  # K4 plot (seasonal)
  lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
         mean(place[9:10]),mean(place[11:12]))
  labelx<-c("11","1","3","5","7","9")
  
  top <- 3
  bottom <- -3
  
  plot(rep(date2,2),scale(K4[i,]),
       type="p",las=1,tcl=-0.2,xlim=c(0,334),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  lines(c(1:334), scale(y.K4[c(305:365,1:273)]), col=rgb(0.9,0,0))
  
  axis(side=1,at=place,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  axis(side=2,at=c(-3,0,3),tcl=-0.1,las=2,pos=0)
  arrows(0,bottom,0,top,length=0)	
  arrows(0,top,334,top,length=0)
  arrows(334,bottom,334,top,length=0)	
  
  
  # K4 plot (diel)
  tp<-c(12,24,48,63)
  lit<-c(mean(tp[1:2]),mean(tp[2:3]),mean(tp[3:4]))
  labelx<-c("13","14","15")
  
  plot(scale(K4.diel[i,]),type="n",las=1,tcl=-0.2,xlim=c(12,63),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.5)
  rect(18.17,bottom,29.7,top,col="gray85",border="transparent")
  rect(42.15,bottom,53.72,top,col="gray85",border="transparent")
  par(new=T)
  plot(time2,scale(K4.diel[i,]),
       type="p",las=1,tcl=-0.2,xlim=c(12,63),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  lines(time3, scale(y.K4.diel), col=rgb(0.9,0,0))
  
  axis(side=1,at=tp,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  axis(side=2,at=c(-3,0,3),tcl=-0.1,las=2,pos=12)
  arrows(12,bottom,12,top,length=0)	
  arrows(12,top,63,top,length=0)
  arrows(63,bottom,63,top,length=0)	
  
  set.panel()
  dev.off()	
}
