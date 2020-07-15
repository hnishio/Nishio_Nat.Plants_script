
library(edgeR)
library(VennDiagram)
library(season)
library(fields)



##### Extraction of genes with seasonal difference #####
# RNA
load("20181120_bowtie2_rpm.Rdata")
RNA.mat <- rawcnt[8410:nrow(rawcnt),1:490]
nrow(RNA.mat)

attribute <- read.table("140519_SampleAttribute.txt", header=T,sep="\t")
attribute.date <- attribute[1:490,2:5]
times <- as.POSIXct(paste(attribute.date[,'year'], 
                          attribute.date[,'month'], 
                          attribute.date[,'day'],
                          attribute.date[,'hour'],
                          sep='/'), 
                    format='%Y/%m/%d/%H')
date <- as.Date(times)
colnames(RNA.mat) <- as.character(date)

Ahal_ID <- rep(NA,length=nrow(RNA.mat))
list <- strsplit(row.names(RNA.mat), ".t")
for(i in 1:length(Ahal_ID)){
  Ahal_ID[i] <- list[[i]][1]
}
row.names(RNA.mat) <- Ahal_ID

head(RNA.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(colnames(RNA.mat))

# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = RNA.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 3
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalRNA_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:length(unique(colnames(RNA.mat))))
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(RNA.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalRNA_CPMover1_edgeR_result.csv", row.names=F)



##### Cosinor model for RNA expression #####

##### Preparation for cosinor #####
load("20181120_bowtie2_rpm.Rdata")
RNA <- log2(rpkm[8410:nrow(rpkm),1:490]+1)
nrow(RNA)

attribute <- read.table("140519_SampleAttribute.txt", header=T,sep="\t")
attribute.date <- attribute[1:490,2:5]
times <- as.POSIXct(paste(attribute.date[,'year'], 
                          attribute.date[,'month'], 
                          attribute.date[,'day'],
                          attribute.date[,'hour'],
                          sep='/'), 
                    format='%Y/%m/%d/%H')
date <- as.Date(times)


RNA.mat <- matrix(nrow=nrow(RNA),ncol=9)

for(i in 1:nrow(RNA)){
  RNA.data <- data.frame(date=date,RNA=RNA[i,])
  
  res.RNA <- cosinor(RNA~1,date="date",data=RNA.data)
  
  phase <- gsub("Month = ", "2012-", summary(res.RNA)$phase)
  phase <- gsub("月 , day = ", "-", phase)
  peak <- as.numeric(as.Date(phase) - as.Date("2012-01-01") + 1)
  
  lphase <- gsub("Month = ", "2012-", summary(res.RNA)$lphase)
  lphase <- gsub("月 , day = ", "-", lphase)
  trough <- as.numeric(as.Date(lphase) - as.Date("2012-01-01") + 1)
  
  RNA.mat[i,] <- c(summary(res.RNA)$amp*2,peak,trough,
                   summary(res.RNA$glm)$coefficients[,"Estimate"],
                   summary(res.RNA$glm)$coefficients[,"Pr(>|t|)"])
}

colnames(RNA.mat) <- c("amp.RNA","peak.RNA","trough.RNA",
                       "coef.inter.RNA","coef.cosw.RNA","coef.sinw.RNA",
                       "pval.inter.RNA","pval.cosw.RNA","pval.sinw.RNA")


Ahal_ID <- rep(NA,length=nrow(RNA))
list <- strsplit(row.names(RNA), ".t")
for(i in 1:length(Ahal_ID)){
  Ahal_ID[i] <- list[[i]][1]
}

RNA.mat <- cbind(Ahal_ID,RNA.mat)
write.csv(RNA.mat, file="RNA_cosinor.csv", quote=F, row.names=F)


# Extraction of expressed genes
RNA.mat <- read.csv("RNA_cosinor.csv",header=T,sep=",")
RNA.CPMover1 <- read.csv("seasonalRNA_CPMover1.csv",header=T,sep=",")
RNA <- merge(RNA.mat,RNA.CPMover1,by="Ahal_ID",all=F,sort=F)
write.csv(RNA[,1:ncol(RNA.mat)], file="RNA_cosinor_CPMover1.csv", quote=F, row.names=F)



##### Comparison between gene list of seasonal ampover3, edgeR, and cosinor #####

# RNA
RNA.over1 <- read.csv("RNA_stat_maxRNAover2_ampRNAover3.csv",header=T,sep=",")
RNA.edgeR <- read.csv("seasonalRNA_CPMover1_edgeR_result.csv",header=T,sep=",")
RNA.cosinor <- read.csv("RNA_cosinor_CPMover1.csv",header=T,sep=",")
nrow(RNA.over1)

length(which(is.na(RNA.edgeR$FDR)))
RNA.edgeRfdr0.05 <- RNA.edgeR[RNA.edgeR$FDR<0.05,]
nrow(RNA.edgeRfdr0.05)

length(which(is.na(RNA.cosinor$pval.cosw.RNA)))
length(which(is.na(RNA.cosinor$pval.sinw.RNA)))
RNA.cosinor0.05 <- RNA.cosinor[(RNA.cosinor$pval.cosw.RNA<0.025 | RNA.cosinor$pval.sinw.RNA<0.025),]
nrow(RNA.cosinor0.05)

#n12
RNA.over1.edgeRfdr0.05 <- merge(RNA.over1, RNA.edgeRfdr0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(RNA.over1.edgeRfdr0.05)

#n123
RNA.over1.edgeRfdr0.05.cosinor0.05 <- merge(RNA.over1.edgeRfdr0.05, RNA.cosinor0.05, by=c("Ahal_ID"),all=F,sort=F)
nrow(RNA.over1.edgeRfdr0.05.cosinor0.05)

write.csv(RNA.over1.edgeRfdr0.05.cosinor0.05[,1:ncol(RNA.over1)], 
          file="RNA_stat_maxRNAover2_ampRNAover3_edgeR_cosinor.csv", row.names=F)
