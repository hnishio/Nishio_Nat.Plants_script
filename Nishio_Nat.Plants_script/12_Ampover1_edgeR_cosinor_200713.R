
library(VennDiagram)



##### Comparison between gene list of seasonal ampover1, edgeR, and cosinor #####

# K4
K4.over1 <- read.csv("K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
K4.edgeR <- read.csv("seasonalK4_CPMover1_edgeR_result.csv",header=T,sep=",")
K4.cosinor <- read.csv("K4K27_seasonal_cosinor_K4_CPMover1.csv",header=T,sep=",")
nrow(K4.over1)

length(which(is.na(K4.edgeR$FDR)))
K4.edgeRfdr0.05 <- K4.edgeR[K4.edgeR$FDR<0.05,]
nrow(K4.edgeRfdr0.05)

length(which(is.na(K4.cosinor$pval.cosw.K4)))
length(which(is.na(K4.cosinor$pval.sinw.K4)))
K4.cosinor0.05 <- K4.cosinor[(K4.cosinor$pval.cosw.K4<0.025 | K4.cosinor$pval.sinw.K4<0.025),]
nrow(K4.cosinor0.05)

#n12
K4.over1.edgeRfdr0.05 <- merge(K4.over1, K4.edgeRfdr0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.over1.edgeRfdr0.05)

#n23
K4.edgeRfdr0.05.cosinor0.05 <- merge(K4.edgeRfdr0.05, K4.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.edgeRfdr0.05.cosinor0.05)

#n13
K4.over1.cosinor0.05 <- merge(K4.over1, K4.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.over1.cosinor0.05)

#n123
K4.over1.edgeRfdr0.05.cosinor0.05 <- merge(K4.over1.edgeRfdr0.05, K4.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.over1.edgeRfdr0.05.cosinor0.05)

colnames(K4.over1.edgeRfdr0.05.cosinor0.05) <- gsub(".x", "", colnames(K4.over1.edgeRfdr0.05.cosinor0.05))
colnames(K4.over1.edgeRfdr0.05.cosinor0.05) <- gsub("mK", "maxK", colnames(K4.over1.edgeRfdr0.05.cosinor0.05))
write.csv(K4.over1.edgeRfdr0.05.cosinor0.05[,1:ncol(K4.over1)], 
          file="K4K27_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv", row.names=F)


# Venn diagram
dir.create("../figs/edgeR")

pdf("../figs/edgeR/venn_K4_amp_FDR_cosinor.pdf",height=2.5,width=2.5)
draw.triple.venn(
  area1=nrow(K4.over1), area2=nrow(K4.edgeRfdr0.05), area3=nrow(K4.cosinor0.05),
  n12=nrow(K4.over1.edgeRfdr0.05), n23=nrow(K4.edgeRfdr0.05.cosinor0.05),
  n13=nrow(K4.over1.cosinor0.05),n123=nrow(K4.over1.edgeRfdr0.05.cosinor0.05),
  cex=1, category=c("K4 over1","K4 FDR < 0.05","K4 cosinor"), cat.cex=c(1,1,1), 
  fontfamily="sans", col=rep("transparent",3),
  fill=c(rgb(0.5,0,0),rgb(1,0,0),"orange"), alpha=rep(0.5,3),
  cat.pos=c(330,30,180), cat.dist=rep(0.05,3), cat.fontfamily="sans",
  cat.col=c(rgb(0.5,0,0),rgb(1,0,0),"orange"), margin=0.1
)
dev.off()


# K27
K27.over1 <- read.csv("K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
K27.edgeR <- read.csv("seasonalK27_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.cosinor <- read.csv("K4K27_seasonal_cosinor_K27_CPMover1.csv",header=T,sep=",")
nrow(K27.over1)

length(which(is.na(K27.edgeR$FDR)))
K27.edgeRfdr0.05 <- K27.edgeR[K27.edgeR$FDR<0.05,]
nrow(K27.edgeRfdr0.05)

length(which(is.na(K27.cosinor$pval.cosw.K27)))
length(which(is.na(K27.cosinor$pval.sinw.K27)))
K27.cosinor0.05 <- K27.cosinor[(K27.cosinor$pval.cosw.K27<0.025 | K27.cosinor$pval.sinw.K27<0.025),]
nrow(K27.cosinor0.05)

#n12
K27.over1.edgeRfdr0.05 <- merge(K27.over1, K27.edgeRfdr0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.over1.edgeRfdr0.05)

#n23
K27.edgeRfdr0.05.cosinor0.05 <- merge(K27.edgeRfdr0.05, K27.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.edgeRfdr0.05.cosinor0.05)

#n13
K27.over1.cosinor0.05 <- merge(K27.over1, K27.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.over1.cosinor0.05)

#n123
K27.over1.edgeRfdr0.05.cosinor0.05 <- merge(K27.over1.edgeRfdr0.05, K27.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.over1.edgeRfdr0.05.cosinor0.05)

colnames(K27.over1.edgeRfdr0.05.cosinor0.05) <- gsub(".x", "", colnames(K27.over1.edgeRfdr0.05.cosinor0.05))
colnames(K27.over1.edgeRfdr0.05.cosinor0.05) <- gsub("mK", "maxK", colnames(K27.over1.edgeRfdr0.05.cosinor0.05))
write.csv(K27.over1.edgeRfdr0.05.cosinor0.05[,1:ncol(K27.over1)], 
          file="K4K27_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv", row.names=F)


# Venn diagram
pdf("../figs/edgeR/venn_K27_amp_FDR_cosinor.pdf",height=2.5,width=2.5)
draw.triple.venn(
  area1=nrow(K27.over1), area2=nrow(K27.edgeRfdr0.05), area3=nrow(K27.cosinor0.05),
  n12=nrow(K27.over1.edgeRfdr0.05), n23=nrow(K27.edgeRfdr0.05.cosinor0.05),
  n13=nrow(K27.over1.cosinor0.05),n123=nrow(K27.over1.edgeRfdr0.05.cosinor0.05),
  cex=1, category=c("K27 over1","K27 FDR < 0.05","K27 cosinor"), cat.cex=c(1,1,1), 
  fontfamily="sans", col=rep("transparent",3),
  fill=c(rgb(0,0,0.3),rgb(0,0,1),"skyblue3"), alpha=rep(0.5,3),
  cat.pos=c(330,30,180), cat.dist=rep(0.05,3), cat.fontfamily="sans",
  cat.col=c(rgb(0,0,0.3),rgb(0,0,1),"skyblue3"), margin=0.1
)
dev.off()



##### Comparison between gene list of diel ampover1, edgeR, and cosinor #####

# K4
K4.over1 <- read.csv("K4K27OMO48h_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")
K4.edgeR <- read.csv("dielK4_CPMover1_edgeR_result.csv",header=T,sep=",")
K4.cosinor <- read.csv("K4K27_diel_cosinor_K4_CPMover1.csv",header=T,sep=",")
nrow(K4.over1)

length(which(is.na(K4.edgeR$FDR)))
K4.edgeRfdr0.05 <- K4.edgeR[K4.edgeR$FDR<0.05,]
nrow(K4.edgeRfdr0.05)

length(which(is.na(K4.cosinor$pval.cosw.K4)))
length(which(is.na(K4.cosinor$pval.sinw.K4)))
K4.cosinor0.05 <- K4.cosinor[(K4.cosinor$pval.cosw.K4<0.025 | K4.cosinor$pval.sinw.K4<0.025),]
nrow(K4.cosinor0.05)

#n12
K4.over1.edgeRfdr0.05 <- merge(K4.over1, K4.edgeRfdr0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.over1.edgeRfdr0.05)

#n23
K4.edgeRfdr0.05.cosinor0.05 <- merge(K4.edgeRfdr0.05, K4.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.edgeRfdr0.05.cosinor0.05)

#n13
K4.over1.cosinor0.05 <- merge(K4.over1, K4.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.over1.cosinor0.05)

#n123
K4.over1.edgeRfdr0.05.cosinor0.05 <- merge(K4.over1.edgeRfdr0.05, K4.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K4.over1.edgeRfdr0.05.cosinor0.05)

colnames(K4.over1.edgeRfdr0.05.cosinor0.05) <- gsub(".x", "", colnames(K4.over1.edgeRfdr0.05.cosinor0.05))
colnames(K4.over1.edgeRfdr0.05.cosinor0.05) <- gsub("mK", "maxK", colnames(K4.over1.edgeRfdr0.05.cosinor0.05))
write.csv(K4.over1.edgeRfdr0.05.cosinor0.05[,1:ncol(K4.over1)], 
          file="K4K27OMO48h_data_stat_maxK4over2_ampK4over1_edgeR_cosinor.csv", row.names=F)


# Venn diagram
pdf("../figs/edgeR/venn_dielK4_amp_FDR_cosinor.pdf",height=2.5,width=2.5)
draw.triple.venn(
  area1=nrow(K4.over1), area2=nrow(K4.edgeRfdr0.05), area3=nrow(K4.cosinor0.05),
  n12=nrow(K4.over1.edgeRfdr0.05), n23=nrow(K4.edgeRfdr0.05.cosinor0.05),
  n13=nrow(K4.over1.cosinor0.05),n123=nrow(K4.over1.edgeRfdr0.05.cosinor0.05),
  cex=1, category=c("K4 over1","K4 FDR < 0.05","K4 cosinor"), cat.cex=c(1,1,1), 
  fontfamily="sans", col=rep("transparent",3),
  fill=c(rgb(0.5,0,0),rgb(1,0,0),"orange"), alpha=rep(0.5,3),
  cat.pos=c(330,30,180), cat.dist=rep(0.05,3), cat.fontfamily="sans",
  cat.col=c(rgb(0.5,0,0),rgb(1,0,0),"orange"), margin=0.1
)
dev.off()


# K27
K27.over1 <- read.csv("K4K27OMO48h_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")
K27.edgeR <- read.csv("dielK27_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.cosinor <- read.csv("K4K27_diel_cosinor_K27_CPMover1.csv",header=T,sep=",")
nrow(K27.over1)

length(which(is.na(K27.edgeR$FDR)))
K27.edgeRfdr0.05 <- K27.edgeR[K27.edgeR$FDR<0.05,]
nrow(K27.edgeRfdr0.05)

length(which(is.na(K27.cosinor$pval.cosw.K27)))
length(which(is.na(K27.cosinor$pval.sinw.K27)))
K27.cosinor0.05 <- K27.cosinor[(K27.cosinor$pval.cosw.K27<0.025 | K27.cosinor$pval.sinw.K27<0.025),]
nrow(K27.cosinor0.05)

#n12
K27.over1.edgeRfdr0.05 <- merge(K27.over1, K27.edgeRfdr0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.over1.edgeRfdr0.05)

#n23
K27.edgeRfdr0.05.cosinor0.05 <- merge(K27.edgeRfdr0.05, K27.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.edgeRfdr0.05.cosinor0.05)

#n13
K27.over1.cosinor0.05 <- merge(K27.over1, K27.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.over1.cosinor0.05)

#n123
K27.over1.edgeRfdr0.05.cosinor0.05 <- merge(K27.over1.edgeRfdr0.05, K27.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(K27.over1.edgeRfdr0.05.cosinor0.05)

colnames(K27.over1.edgeRfdr0.05.cosinor0.05) <- gsub(".x", "", colnames(K27.over1.edgeRfdr0.05.cosinor0.05))
colnames(K27.over1.edgeRfdr0.05.cosinor0.05) <- gsub("mK", "maxK", colnames(K27.over1.edgeRfdr0.05.cosinor0.05))
write.csv(K27.over1.edgeRfdr0.05.cosinor0.05[,1:ncol(K27.over1)], 
          file="K4K27OMO48h_data_stat_maxK27over2_ampK27over1_edgeR_cosinor.csv", row.names=F)


# Venn diagram
pdf("../figs/edgeR/venn_dielK27_amp_FDR_cosinor.pdf",height=2.5,width=2.5)
draw.triple.venn(
  area1=nrow(K27.over1), area2=nrow(K27.edgeRfdr0.05), area3=nrow(K27.cosinor0.05),
  n12=nrow(K27.over1.edgeRfdr0.05), n23=nrow(K27.edgeRfdr0.05.cosinor0.05),
  n13=nrow(K27.over1.cosinor0.05),n123=nrow(K27.over1.edgeRfdr0.05.cosinor0.05),
  cex=1, category=c("K27 over1","K27 FDR < 0.05","K27 cosinor"), cat.cex=c(1,1,1), 
  fontfamily="sans", col=rep("transparent",3),
  fill=c(rgb(0,0,0.3),rgb(0,0,1),"skyblue3"), alpha=rep(0.5,3),
  cat.pos=c(330,30,180), cat.dist=rep(0.05,3), cat.fontfamily="sans",
  cat.col=c(rgb(0,0,0.3),rgb(0,0,1),"skyblue3"), margin=0.1
)
dev.off()
