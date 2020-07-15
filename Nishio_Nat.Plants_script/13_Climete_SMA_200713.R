
library(TTR)
library(fields)



##### Calculation of Simple Moving Average #####
setwd("..")

climate <- read.csv("data/Nishiwaki_data.csv",header=T,sep=",")
plot(climate$day_length, type="l")

# SMA at sampling date
climate.sma <- matrix(ncol=2*4,nrow=12)
sma.period <- c(7,30)
samplepoint <- c(311,339,374,402,430,458,486,514,549,577,605,633)
samplepoint2 <- samplepoint-1
length(samplepoint2)

# Calculation of SMA
for(i in 1:2){
  mvtemp<-SMA(climate$mean_temp,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i] <- mvtemp[samplepoint2[j]]
  }
  
  mvprep<-SMA(climate$precipitation,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i+2] <- mvprep[samplepoint2[j]]
  }
  
  mvday<-SMA(climate$day_length,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i+4] <- mvday[samplepoint2[j]]
  }
  
  mvsun<-SMA(climate$sunlight_hours,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i+6] <- mvsun[samplepoint2[j]]
  }
}

colnames(climate.sma) <- c(paste0("SMA",sma.period,c("_temp")),
                           paste0("SMA",sma.period,c("_ppt")),
                           paste0("SMA",sma.period,c("_day")),
                           paste0("SMA",sma.period,c("_sun"))
)
date <- c("Nov.6","Dec.4","Jan.8","Feb.5","Mar.5","Apr.2",
          "Apr.30","May.28","Jul.2","Jul.30","Aug.27","Sep.24")
climate.sma2 <- cbind(date,climate.sma)

write.csv(climate.sma2,"data/SMA_climate.csv",quote=F,row.names=F)



####### SMA plot ################

# SMA for visualisation (121001-130930)
climate.sma <- matrix(ncol=2*4,nrow=365)
sma.period <- c(7,30)
samplepoint <- c(275:639)
samplepoint2 <- samplepoint-1
length(samplepoint2)


# Calculation of SMA
for(i in 1:2){
  mvtemp<-SMA(climate$mean_temp,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i] <- mvtemp[samplepoint2[j]]
  }
  
  mvprep<-SMA(climate$precipitation,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i+2] <- mvprep[samplepoint2[j]]
  }
  
  mvday<-SMA(climate$day_length,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i+4] <- mvday[samplepoint2[j]]
  }
  
  mvsun<-SMA(climate$sunlight_hours,sma.period[i])
  for(j in 1:length(samplepoint2)){
    climate.sma[j,i+6] <- mvsun[samplepoint2[j]]
  }
}

colnames(climate.sma) <- c(paste0("SMA",sma.period,c("_temp")),
                           paste0("SMA",sma.period,c("_ppt")),
                           paste0("SMA",sma.period,c("_day")),
                           paste0("SMA",sma.period,c("_sun"))
)

date <- as.character(climate$date[275:639])
climate.sma2 <- cbind(date,climate.sma)

write.csv(climate.sma2,"data/SMA_121001-130930.csv",quote=F,row.names=F)


# Graph
x<-read.csv("data/SMA_121001-130930.csv",header=T,sep=",")
tempcol <- c("black","red")

dir.create("figs/climate_correlation/")
pdf("figs/climate_correlation/SMA_climate.pdf", height=1, width=5*4/3)
set.panel(1,4)
par(oma=c(0,0,0,0))
par(mar=c(1.2,1.7,0.8,0.2))
par(mgp=c(0,0.1,0))
par(ps=6)
par(xpd=T)
par(cex=1)

# temperature
plot(1:365,x$SMA7_temp,type="l",pch=0,col=tempcol[1],ann=F,axes=F,
     xlim=c(0,365),ylim=c(-5,35),cex=0.7)
par(new=T)
plot(1:365,x$SMA30_temp,type="l",pch=15,col=tempcol[2],ann=F,axes=F,
     xlim=c(0,365),ylim=c(-5,35),cex=0.7)

day<-c(0,31,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),
       mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),
       mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]))
labelx<-c("10","11","12","1","2","3","4","5","6","7","8","9")
axis(side=1,at=place,label=F,tcl=-0.1,pos=-5)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(-5,35,10),tcl=0.1,las=2,pos=0)
mtext("SMA of",side=2,line=0.8)
mtext("temperature (°C)",side=2,line=0.38)
arrows(0,35,365,35,length=0)
arrows(365,-5,365,35,length=0)

arrows(0,-5-40*0.225,92,-5-40*0.225,length=0.03,code=3,lwd=0.75)
arrows(92,-5-40*0.225,365,-5-40*0.225,length=0.03,code=3,lwd=0.75)
mtext("2012",side=1,at=46,line=0.2)
mtext("2013",side=1,at=229,line=0.2)

legend(-5, 35+40*0.5, "1 week", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[1], pt.cex=0.7)
legend(155, 35+40*0.5, "1 month", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[2], pt.cex=0.7)


# precipitation
plot(1:365,x$SMA7_ppt,type="l",pch=0,col=tempcol[1],ann=F,axes=F,
     xlim=c(0,365),ylim=c(0,40),cex=0.7)
par(new=T)
plot(1:365,x$SMA30_ppt,type="l",pch=15,col=tempcol[2],ann=F,axes=F,
     xlim=c(0,365),ylim=c(0,40),cex=0.7)

day<-c(0,31,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),
       mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),
       mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]))
labelx<-c("10","11","12","1","2","3","4","5","6","7","8","9")
axis(side=1,at=place,label=F,tcl=-0.1,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,40,10),tcl=0.1,las=2,pos=0)
mtext("SMA of",side=2,line=0.8)
mtext("precipitation (mm)",side=2,line=0.38)
arrows(0,40,365,40,length=0)
arrows(365,0,365,40,length=0)

arrows(0,0-40*0.225,92,0-40*0.225,length=0.03,code=3,lwd=0.75)
arrows(92,0-40*0.225,365,0-40*0.225,length=0.03,code=3,lwd=0.75)
mtext("2012",side=1,at=46,line=0.2)
mtext("2013",side=1,at=229,line=0.2)

legend(-5, 40+40*0.5, "1 week", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[1], pt.cex=0.7)
legend(155, 40+40*0.5, "1 month", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[2], pt.cex=0.7)


# day length
plot(1:365,x$SMA7_day,type="l",pch=0,col=tempcol[1],ann=F,axes=F,
     xlim=c(0,365),ylim=c(8,16),cex=0.7)
par(new=T)
plot(1:365,x$SMA30_day,type="l",pch=15,col=tempcol[2],ann=F,axes=F,
     xlim=c(0,365),ylim=c(8,16),cex=0.7)

day<-c(0,31,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),
       mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),
       mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]))
labelx<-c("10","11","12","1","2","3","4","5","6","7","8","9")
axis(side=1,at=place,label=F,tcl=-0.1,pos=8)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(8,16,4),tcl=0.1,las=2,pos=0)
mtext("SMA of",side=2,line=0.8)
mtext("day length (hours)",side=2,line=0.38)
arrows(0,16,365,16,length=0)
arrows(365,8,365,16,length=0)

arrows(0,8-8*0.225,92,8-8*0.225,length=0.03,code=3,lwd=0.75)
arrows(92,8-8*0.225,365,8-8*0.225,length=0.03,code=3,lwd=0.75)
mtext("2012",side=1,at=46,line=0.2)
mtext("2013",side=1,at=229,line=0.2)

legend(-5, 16+8*0.5, "1 week", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[1], pt.cex=0.7)
legend(155, 16+8*0.5, "1 month", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[2], pt.cex=0.7)


# sunlight hours
plot(1:365,x$SMA7_sun,type="l",pch=0,col=tempcol[1],ann=F,axes=F,
     xlim=c(0,365),ylim=c(0,10),cex=0.7)
par(new=T)
plot(1:365,x$SMA30_sun,type="l",pch=15,col=tempcol[2],ann=F,axes=F,
     xlim=c(0,365),ylim=c(0,10),cex=0.7)

day<-c(0,31,30,31,31,28,31,30,31,30,31,31,30)
place<-c(day[1],sum(day[1:2]),sum(day[1:3]),sum(day[1:4]),sum(day[1:5]),
         sum(day[1:6]),sum(day[1:7]),sum(day[1:8]),sum(day[1:9]),sum(day[1:10]),
         sum(day[1:11]),sum(day[1:12]),sum(day[1:13]))
lit<-c(mean(place[1:2]),mean(place[2:3]),mean(place[3:4]),mean(place[4:5]),
       mean(place[5:6]),mean(place[6:7]),mean(place[7:8]),mean(place[8:9]),
       mean(place[9:10]),mean(place[10:11]),mean(place[11:12]),mean(place[12:13]))
labelx<-c("10","11","12","1","2","3","4","5","6","7","8","9")
axis(side=1,at=place,label=F,tcl=-0.1,pos=0)
mtext(side=1,at=lit,labelx,line=-0.45)
axis(side=2,at=seq(0,10,5),tcl=0.1,las=2,pos=0)
mtext("SMA of",side=2,line=0.8)
mtext("sunlight hours",side=2,line=0.38)
arrows(0,10,365,10,length=0)
arrows(365,0,365,10,length=0)

arrows(0,0-10*0.225,92,0-10*0.225,length=0.03,code=3,lwd=0.75)
arrows(92,0-10*0.225,365,0-10*0.225,length=0.03,code=3,lwd=0.75)
mtext("2012",side=1,at=46,line=0.2)
mtext("2013",side=1,at=229,line=0.2)

legend(-5, 10+10*0.5, "1 week", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[1], pt.cex=0.7)
legend(155, 10+10*0.5, "1 month", lwd=1, x.intersp=0.1, seg.len=0.7, bty="n", 
       col=tempcol[2], pt.cex=0.7)

dev.off()
