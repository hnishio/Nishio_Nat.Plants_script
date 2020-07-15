
library(fields)



##### Graphs of seasonal and diel change #####
setwd("data")

###### Preparation for spline #####

# Load the seasonal data

K4.K27 <- read.csv("K4K27_data_stat.csv",
                   header=T,sep=",")
K4.K27 <- K4.K27[order(K4.K27$Ahal_ID),]

K4.rep1 <- cbind(K4.K27[,11:22],K4.K27[,11:22],K4.K27[,11:22])
K27.rep1 <- cbind(K4.K27[,26:37],K4.K27[,26:37],K4.K27[,26:37])
K4.rep2 <- cbind(K4.K27[,41:52],K4.K27[,41:52],K4.K27[,41:52])
K27.rep2 <- cbind(K4.K27[,56:67],K4.K27[,56:67],K4.K27[,56:67])
K4 <- cbind(K4.rep1,K4.rep2)
K27 <- cbind(K27.rep1,K27.rep2)
K4 <- as.matrix(K4)
K27 <- as.matrix(K27)
date <- c(6, 34, 69, 97, 125, 153, 181, 209, 244, 272, 300, 328)
date <- c(date,365+date,730+date)
date <- c(date,date)
title <- K4.K27[,c(1:6,37)]

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

K4.rep1 <- cbind(K4.K27.rep1[,11:18],K4.K27.rep1[,11:18],K4.K27.rep1[,11:18])
K27.rep1 <- cbind(K4.K27.rep1[,22:29],K4.K27.rep1[,22:29],K4.K27.rep1[,22:29])
K4.rep2 <- cbind(K4.K27.rep2[,11:18],K4.K27.rep2[,11:18],K4.K27.rep2[,11:18])
K27.rep2 <- cbind(K4.K27.rep2[,22:29],K4.K27.rep2[,22:29],K4.K27.rep2[,22:29])
K4.rep3 <- cbind(K4.K27.rep3[,11:18],K4.K27.rep3[,11:18],K4.K27.rep3[,11:18])
K27.rep3 <- cbind(K4.K27.rep3[,22:29],K4.K27.rep3[,22:29],K4.K27.rep3[,22:29])
K4.rep4 <- cbind(K4.K27.rep4[,11:18],K4.K27.rep4[,11:18],K4.K27.rep4[,11:18])
K27.rep4 <- cbind(K4.K27.rep4[,22:29],K4.K27.rep4[,22:29],K4.K27.rep4[,22:29])

K4.diel <- cbind(K4.rep1,K4.rep2,K4.rep3,K4.rep4)
K27.diel <- cbind(K27.rep1,K27.rep2,K27.rep3,K27.rep4)
K4.diel <- as.matrix(K4.diel)
K27.diel <- as.matrix(K27.diel)
time <- c(18, 24, 30, 36, 42, 48, 54, 60)
time <- c(time,48+time,96+time)
time <- c(time,time,time,time)



##### Seasonal change #####
dir.create("../figs/seasonal_diel_change")
day<-c(365,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]))


for(i in c(1263,1268,1821,2309,7214,8236,9640,9784,9826,10555,10811,13729,14931,15060,16823,17199,17215,17335,18162,18589,19190,20174,20212,20684,27591,28191)){
  pdf(paste("../figs/seasonal_diel_change/K27K4seasonal_diel_", 
            title$Ahal_ID[i], "_", title$Symbol[i], ".pdf", sep=""),
      height=0.45,width=2.95) 
  set.panel(1,4,relax=T)
  par(oma=c(0,0,0,0))
  par(mar=c(0.6,0.3,0.2,0.1))
  par(mgp=c(3,0.13,0.15))
  par(ps=6)
  par(cex=1)
  

# Seasonal change
# K27 plot
  lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
       mean(place[9:10]),mean(place[11:12]))
  labelx<-c("11","1","3","5","7","9")
  
  spK27 <- smooth.spline(date,K27[i,],spar=0.3)
  length <- 1000
  x <- seq(365,699,length=length)
  predK27 <- predict(spK27,x)
  
  if((max(c(K27[i,],K27.diel[i,])) + max(K27[i,])*0.1) >= 1){
  top <- max(c(K27[i,],K27.diel[i,])) + max(K27[i,])*0.1
  }else{top <- 1}
  if((min(c(K27[i,],K27.diel[i,])) - max(c(K27[i,],K27.diel[i,]))*0.1) > (max(c(K27[i,],K27.diel[i,]))-min(c(K27[i,],K27.diel[i,])))/10){
    bottom <- min(c(K27[i,],K27.diel[i,])) - max(c(K27[i,],K27.diel[i,]))*0.1
    if(bottom < 1){bottom <- 0}
  }else{bottom <- 0}
  
  plot(date[c(13:24,49:60)],K27[i,c(13:24,49:60)],
       type="p",las=1,tcl=-0.2,xlim=c(365,699),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  par(new=T)
  plot(predK27, type="l", col=rgb(0,0,0.7),xlim=c(365,699),ylim=c(bottom,top),
       xlab="",ylab="",ann=F,axes=F)
  
  axis(side=1,at=place,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  
  
  if(bottom==0){
    if((floor(top)+floor(bottom))%%2==0){
      axis(side=2,at=c(0,
                       mean(c(floor(bottom),floor(top))),
                       floor(top)),
           tcl=-0.1,las=2,pos=365)
    }else{if(floor(top)==1){axis(side=2,at=c(0,1),tcl=-0.1,las=2,pos=365)}else{
      axis(side=2,at=c(0, (floor(top)-1)/2, floor(top)-1),
           tcl=-0.1,las=2,pos=365)}
    }}else{
      
      #if(floor(top)%%2==0){
      if((floor(top)+floor(bottom))%%2==0){
        axis(side=2,at=c(floor(bottom),
                         mean(c(floor(bottom),floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=365)
      }else{
        axis(side=2,at=c(floor(bottom)+1,
                         mean(c(floor(bottom)+1,floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=365)
      }
    }
  
  arrows(365,bottom,365,top,length=0)	
  arrows(365,top,699,top,length=0)
  arrows(699,bottom,699,top,length=0)	


# diel change
# K27 plot
  tp<-c(60,72,96,114)
  lit<-c(mean(tp[1:2]),mean(tp[2:3]),mean(tp[3:4]))
  labelx<-c("13","14","15")
  
  spK27 <- smooth.spline(time,K27.diel[i,],spar=0.3)
  length <- 1000
  x <- seq(60,114,length=length)
  predK27 <- predict(spK27,x)
  
  if((max(c(K27[i,],K27.diel[i,])) + max(K27[i,])*0.1) >= 1){
  top <- max(c(K27[i,],K27.diel[i,])) + max(K27[i,])*0.1
  }else{top <- 1}
  if((min(c(K27[i,],K27.diel[i,])) - max(c(K27[i,],K27.diel[i,]))*0.1) > (max(c(K27[i,],K27.diel[i,]))-min(c(K27[i,],K27.diel[i,])))/10){
    bottom <- min(c(K27[i,],K27.diel[i,])) - max(c(K27[i,],K27.diel[i,]))*0.1
    if(bottom < 1){bottom <- 0}
  }else{bottom <- 0}
  
  plot(predK27,type="n",las=1,tcl=-0.2,xlim=c(60,114),ylim=c(bottom,top),
  	   xlab="",ylab="",axes=F, ann=F, cex=0.5)
  rect(66.15,0,77.68,top,col="gray85",border="transparent")
  rect(90.13,0,101.7,top,col="gray85",border="transparent")
  par(new=T)
  plot(time[c(9:16,33:40,57:64,81:88)],K27.diel[i,c(9:16,33:40,57:64,81:88)],
       type="p",las=1,tcl=-0.2,xlim=c(60,114),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0,0,0.7))
  par(new=T)
  plot(predK27, type="l", col=rgb(0,0,0.7),xlim=c(60,114),ylim=c(bottom,top),
       xlab="",ylab="",ann=F,axes=F)
  
  axis(side=1,at=tp,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  
  
  if(bottom==0){
    if((floor(top)+floor(bottom))%%2==0){
      axis(side=2,at=c(0,
                       mean(c(floor(bottom),floor(top))),
                       floor(top)),
           tcl=-0.1,las=2,pos=60)
    }else{if(floor(top)==1){axis(side=2,at=c(0,1),tcl=-0.1,las=2,pos=60)}else{
      axis(side=2,at=c(0, (floor(top)-1)/2, floor(top)-1),
           tcl=-0.1,las=2,pos=60)}
    }}else{
      
      #if(floor(top)%%2==0){
      if((floor(top)+floor(bottom))%%2==0){
        axis(side=2,at=c(floor(bottom),
                         mean(c(floor(bottom),floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=60)
      }else{if(floor(top)==1){axis(side=2,at=c(0,1),tcl=-0.1,las=2,pos=60)}else{
        axis(side=2,at=c(floor(bottom)+1,
                         mean(c(floor(bottom)+1,floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=60)}
      }
    }
  
  arrows(60,bottom,60,top,length=0)	
  arrows(60,top,114,top,length=0)
  arrows(114,bottom,114,top,length=0)	


# Seasonal change
# K4 plot
  lit<-c(mean(place[1:2]),mean(place[3:4]),mean(place[5:6]),mean(place[7:8]),
       mean(place[9:10]),mean(place[11:12]))
  labelx<-c("11","1","3","5","7","9")
  
  spK4 <- smooth.spline(date,K4[i,],spar=0.3)
  length <- 1000
  x <- seq(365,699,length=length)
  predK4 <- predict(spK4,x)
  
  if((max(c(K4[i,],K4.diel[i,])) + max(K4[i,])*0.1) >= 1){
  top <- max(c(K4[i,],K4.diel[i,])) + max(K4[i,])*0.1
  }else{top <- 1}
  if((min(c(K4[i,],K4.diel[i,])) - max(c(K4[i,],K4.diel[i,]))*0.1) > (max(c(K4[i,],K4.diel[i,]))-min(c(K4[i,],K4.diel[i,])))/10){
    bottom <- min(c(K4[i,],K4.diel[i,])) - max(c(K4[i,],K4.diel[i,]))*0.1
    if(bottom < 1){bottom <- 0}
  }else{bottom <- 0}
  
  plot(date[c(13:24,49:60)],K4[i,c(13:24,49:60)],
       type="p",las=1,tcl=-0.2,xlim=c(365,699),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  par(new=T)
  plot(predK4, type="l", col=rgb(0.9,0,0),xlim=c(365,699),ylim=c(bottom,top),
       xlab="",ylab="",ann=F,axes=F)
  
  axis(side=1,at=place,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  
  
  if(bottom==0){
    if((floor(top)+floor(bottom))%%2==0){
      axis(side=2,at=c(0,
                       mean(c(floor(bottom),floor(top))),
                       floor(top)),
           tcl=-0.1,las=2,pos=365)
    }else{if(floor(top)==1){axis(side=2,at=c(0,1),tcl=-0.1,las=2,pos=365)}else{
      axis(side=2,at=c(0, (floor(top)-1)/2, floor(top)-1),
           tcl=-0.1,las=2,pos=365)}
    }}else{
      
      #if(floor(top)%%2==0){
      if((floor(top)+floor(bottom))%%2==0){
        axis(side=2,at=c(floor(bottom),
                         mean(c(floor(bottom),floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=365)
      }else{
        axis(side=2,at=c(floor(bottom)+1,
                         mean(c(floor(bottom)+1,floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=365)
      }
    }
  
  arrows(365,bottom,365,top,length=0)	
  arrows(365,top,699,top,length=0)
  arrows(699,bottom,699,top,length=0)	


# diel change
# K4 plot
  tp<-c(60,72,96,114)
  lit<-c(mean(tp[1:2]),mean(tp[2:3]),mean(tp[3:4]))
  labelx<-c("13","14","15")
  
  spK4 <- smooth.spline(time,K4.diel[i,],spar=0.3)
  length <- 1000
  x <- seq(60,114,length=length)
  predK4 <- predict(spK4,x)
  
  if((max(c(K4[i,],K4.diel[i,])) + max(K4[i,])*0.1) >= 1){
  top <- max(c(K4[i,],K4.diel[i,])) + max(K4[i,])*0.1
  }else{top <- 1}
  if((min(c(K4[i,],K4.diel[i,])) - max(c(K4[i,],K4.diel[i,]))*0.1) > (max(c(K4[i,],K4.diel[i,]))-min(c(K4[i,],K4.diel[i,])))/10){
    bottom <- min(c(K4[i,],K4.diel[i,])) - max(c(K4[i,],K4.diel[i,]))*0.1
    if(bottom < 1){bottom <- 0}
  }else{bottom <- 0}
  
  plot(predK4,type="n",las=1,tcl=-0.2,xlim=c(60,114),ylim=c(bottom,top),
  	   xlab="",ylab="",axes=F, ann=F, cex=0.5)
  rect(66.15,0,77.68,top,col="gray85",border="transparent")
  rect(90.13,0,101.7,top,col="gray85",border="transparent")
  par(new=T)
  plot(time[c(9:16,33:40,57:64,81:88)],K4.diel[i,c(9:16,33:40,57:64,81:88)],
       type="p",las=1,tcl=-0.2,xlim=c(60,114),ylim=c(bottom,top),
       xlab="",ylab="",axes=F, ann=F, cex=0.3, col=rgb(0.9,0,0))
  par(new=T)
  plot(predK4, type="l", col=rgb(0.9,0,0),xlim=c(60,114),ylim=c(bottom,top),
       xlab="",ylab="",ann=F,axes=F)
  
  axis(side=1,at=tp,label=F,tcl=-0.1,pos=bottom)
  mtext(side=1,at=lit,labelx,line=-0.4)
  
  
  if(bottom==0){
    if((floor(top)+floor(bottom))%%2==0){
      axis(side=2,at=c(0,
                       mean(c(floor(bottom),floor(top))),
                       floor(top)),
           tcl=-0.1,las=2,pos=60)
    }else{if(floor(top)==1){axis(side=2,at=c(0,1),tcl=-0.1,las=2,pos=60)}else{
      axis(side=2,at=c(0, (floor(top)-1)/2, floor(top)-1),
           tcl=-0.1,las=2,pos=60)}
    }}else{
      
      #if(floor(top)%%2==0){
      if((floor(top)+floor(bottom))%%2==0){
        axis(side=2,at=c(floor(bottom),
                         mean(c(floor(bottom),floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=60)
      }else{
        axis(side=2,at=c(floor(bottom)+1,
                         mean(c(floor(bottom)+1,floor(top))),
                         floor(top)),
             tcl=-0.1,las=2,pos=60)
      }
    }
  
  arrows(60,bottom,60,top,length=0)	
  arrows(60,top,114,top,length=0)
  arrows(114,bottom,114,top,length=0)	


  dev.off()	
}

