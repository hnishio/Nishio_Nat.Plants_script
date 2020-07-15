
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)
library(VennDiagram)



##### Extraction of genes with seasonal difference #####
setwd("data")

# H3K4me3
K4.rep1 <- read.csv("K4_rep1.csv",header=T,sep=",")
K4.rep2 <- read.csv("K4_rep2.csv",header=T,sep=",")

colnames(K4.rep1) <- c("K4_chr","K4_start","K4_end", 
                       "K4_rep1_11.6","K4_rep1_12.4","K4_rep1_1.8","K4_rep1_2.5",
                       "K4_rep1_3.5","K4_rep1_4.2","K4_rep1_4.30","K4_rep1_5.28",
                       "K4_rep1_7.2","K4_rep1_7.30","K4_rep1_8.27","K4_rep1_9.24")
colnames(K4.rep2) <- c("K4_chr","K4_start", "K4_end", 
                       "K4_rep2_11.6","K4_rep2_12.4","K4_rep2_1.8","K4_rep2_2.5",
                       "K4_rep2_3.5","K4_rep2_4.2","K4_rep2_4.30","K4_rep2_5.28",
                       "K4_rep2_7.2","K4_rep2_7.30","K4_rep2_8.27","K4_rep2_9.24")
K4.rep1.rep2 <- merge(K4.rep1, K4.rep2, 
                      by=c("K4_chr","K4_start","K4_end"), 
                      all=F, sort=F)

gene <- read.table("TSSdown1000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K4_chr","K4_start","K4_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K4.rep1.rep2.gene <- merge(K4.rep1.rep2, gene, by=c("K4_chr","K4_start","K4_end"), 
                           all=F, sort=F)
K4.rep1.rep2.gene <- K4.rep1.rep2.gene[!duplicated(paste(K4.rep1.rep2.gene$K4_chr,K4.rep1.rep2.gene$K4_start,
                                                         K4.rep1.rep2.gene$K4_end,K4.rep1.rep2.gene$Ahal_ID),sep=","),]
K4.rep1.rep2.gene <- K4.rep1.rep2.gene[,-((ncol(K4.rep1.rep2.gene)-5):ncol(K4.rep1.rep2.gene))]

K4.mat <- K4.rep1.rep2.gene[,4:27]
row.names(K4.mat) <- K4.rep1.rep2.gene[,28]
head(K4.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K4.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK4_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K4.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK4_CPMover1_edgeR_result.csv", row.names=F)


# H3K27me3
K27.rep1 <- read.csv("K27_rep1.csv",header=T,sep=",")
K27.rep2 <- read.csv("K27_rep2.csv",header=T,sep=",")

colnames(K27.rep1) <- c("K27_chr","K27_start","K27_end", 
                       "K27_rep1_11.6","K27_rep1_12.4","K27_rep1_1.8","K27_rep1_2.5",
                       "K27_rep1_3.5","K27_rep1_4.2","K27_rep1_4.30","K27_rep1_5.28",
                       "K27_rep1_7.2","K27_rep1_7.30","K27_rep1_8.27","K27_rep1_9.24")
colnames(K27.rep2) <- c("K27_chr","K27_start", "K27_end", 
                       "K27_rep2_11.6","K27_rep2_12.4","K27_rep2_1.8","K27_rep2_2.5",
                       "K27_rep2_3.5","K27_rep2_4.2","K27_rep2_4.30","K27_rep2_5.28",
                       "K27_rep2_7.2","K27_rep2_7.30","K27_rep2_8.27","K27_rep2_9.24")
K27.rep1.rep2 <- merge(K27.rep1, K27.rep2, 
                      by=c("K27_chr","K27_start","K27_end"), 
                      all=F, sort=F)

gene <- read.table("Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K27_chr","K27_start","K27_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K27.rep1.rep2.gene <- merge(K27.rep1.rep2, gene, by=c("K27_chr","K27_start","K27_end"), 
                           all=F, sort=F)
K27.rep1.rep2.gene <- K27.rep1.rep2.gene[!duplicated(paste(K27.rep1.rep2.gene$K27_chr,K27.rep1.rep2.gene$K27_start,
                                                         K27.rep1.rep2.gene$K27_end,K27.rep1.rep2.gene$Ahal_ID),sep=","),]
K27.rep1.rep2.gene <- K27.rep1.rep2.gene[,-((ncol(K27.rep1.rep2.gene)-5):ncol(K27.rep1.rep2.gene))]

K27.mat <- K27.rep1.rep2.gene[,4:27]
row.names(K27.mat) <- K27.rep1.rep2.gene[,28]
head(K27.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),2))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K27.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="seasonalK27_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:12)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K27.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="seasonalK27_CPMover1_edgeR_result.csv", row.names=F)


# Venn diagram
dir.create("../figs/edgeR")
K4.mat <- read.csv("seasonalK4_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.mat <- read.csv("seasonalK27_CPMover1_edgeR_result.csv",header=T,sep=",")

K4_season_osci <- K4.mat[which(K4.mat$FDR<0.05),]
K27_season_osci <- K27.mat[which(K27.mat$FDR<0.05),]
K4K27_season_osci <- merge(K4_season_osci,K27_season_osci,by=c("Ahal_ID"),all=F,sort=F)


pdf("../figs/edgeR/venn_K4K27_season_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4_season_osci), area2=nrow(K27_season_osci), cross.area=nrow(K4K27_season_osci), 
  cex=1, category=c("K4 seasonal","K27 seasonal"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans", 
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()



##### Extraction of genes with diel difference #####
# H3K4me3
K4.rep1 <- read.csv("K4_rep1_OMO48h.csv", sep=",", header=T)
K4.rep2 <- read.csv("K4_rep2_OMO48h.csv", sep=",", header=T)
K4.rep3 <- read.csv("K4_rep3_OMO48h.csv", sep=",", header=T)
K4.rep4 <- read.csv("K4_rep4_OMO48h.csv", sep=",", header=T)

colnames(K4.rep1) <- c("K4_chr","K4_start","K4_end",
                       "K4_rep1_D1T18","K4_rep1_D1T24","K4_rep1_D2T6","K4_rep1_D2T12",
                       "K4_rep1_D2T18","K4_rep1_D2T24","K4_rep1_D3T6","K4_rep1_D3T12")
colnames(K4.rep2) <- c("K4_chr","K4_start","K4_end",
                       "K4_rep2_D1T18","K4_rep2_D1T24","K4_rep2_D2T6","K4_rep2_D2T12",
                       "K4_rep2_D2T18","K4_rep2_D2T24","K4_rep2_D3T6","K4_rep2_D3T12")
colnames(K4.rep3) <- c("K4_chr","K4_start","K4_end",
                       "K4_rep3_D1T18","K4_rep3_D1T24","K4_rep3_D2T6","K4_rep3_D2T12",
                       "K4_rep3_D2T18","K4_rep3_D2T24","K4_rep3_D3T6","K4_rep3_D3T12")
colnames(K4.rep4) <- c("K4_chr","K4_start","K4_end",
                       "K4_rep4_D1T18","K4_rep4_D1T24","K4_rep4_D2T6","K4_rep4_D2T12",
                       "K4_rep4_D2T18","K4_rep4_D2T24","K4_rep4_D3T6","K4_rep4_D3T12")
K4.rep1.rep2 <- merge(K4.rep1, K4.rep2, 
                      by=c("K4_chr","K4_start","K4_end"), 
                      all=F, sort=F)
K4.rep1.rep2.rep3 <- merge(K4.rep1.rep2, K4.rep3,
                      by=c("K4_chr","K4_start","K4_end"), 
                      all=F, sort=F)
K4.rep1.rep2.rep3.rep4 <- merge(K4.rep1.rep2.rep3, K4.rep4,
                           by=c("K4_chr","K4_start","K4_end"), 
                           all=F, sort=F)

gene <- read.table("TSSdown1000_Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K4_chr","K4_start","K4_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K4.rep1.rep2.rep3.rep4.gene <- merge(K4.rep1.rep2.rep3.rep4, gene, 
                                     by=c("K4_chr","K4_start","K4_end"), 
                                     all=F, sort=F)
K4.rep1.rep2.rep3.rep4.gene <- 
  K4.rep1.rep2.rep3.rep4.gene[!duplicated(paste(K4.rep1.rep2.rep3.rep4.gene$K4_chr,
                                                K4.rep1.rep2.rep3.rep4.gene$K4_start,
                                                K4.rep1.rep2.rep3.rep4.gene$K4_end,
                                                K4.rep1.rep2.rep3.rep4.gene$Ahal_ID),sep=","),]
K4.rep1.rep2.rep3.rep4.gene <- 
  K4.rep1.rep2.rep3.rep4.gene[,-((ncol(K4.rep1.rep2.rep3.rep4.gene)-5):
                                   ncol(K4.rep1.rep2.rep3.rep4.gene))]

K4.mat <- K4.rep1.rep2.rep3.rep4.gene[,4:35]
row.names(K4.mat) <- K4.rep1.rep2.rep3.rep4.gene[,36]
head(K4.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H"),4))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K4.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 4
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="dielK4_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:8)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K4.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="dielK4_CPMover1_edgeR_result.csv", row.names=F)


# H3K27me3
K27.rep1 <- read.csv("K27_rep1_OMO48h.csv", sep=",", header=T)
K27.rep2 <- read.csv("K27_rep2_OMO48h.csv", sep=",", header=T)
K27.rep3 <- read.csv("K27_rep3_OMO48h.csv", sep=",", header=T)
K27.rep4 <- read.csv("K27_rep4_OMO48h.csv", sep=",", header=T)

colnames(K27.rep1) <- c("K27_chr","K27_start","K27_end",
                       "K27_rep1_D1T18","K27_rep1_D1T24","K27_rep1_D2T6","K27_rep1_D2T12",
                       "K27_rep1_D2T18","K27_rep1_D2T24","K27_rep1_D3T6","K27_rep1_D3T12")
colnames(K27.rep2) <- c("K27_chr","K27_start","K27_end",
                       "K27_rep2_D1T18","K27_rep2_D1T24","K27_rep2_D2T6","K27_rep2_D2T12",
                       "K27_rep2_D2T18","K27_rep2_D2T24","K27_rep2_D3T6","K27_rep2_D3T12")
colnames(K27.rep3) <- c("K27_chr","K27_start","K27_end",
                       "K27_rep3_D1T18","K27_rep3_D1T24","K27_rep3_D2T6","K27_rep3_D2T12",
                       "K27_rep3_D2T18","K27_rep3_D2T24","K27_rep3_D3T6","K27_rep3_D3T12")
colnames(K27.rep4) <- c("K27_chr","K27_start","K27_end",
                       "K27_rep4_D1T18","K27_rep4_D1T24","K27_rep4_D2T6","K27_rep4_D2T12",
                       "K27_rep4_D2T18","K27_rep4_D2T24","K27_rep4_D3T6","K27_rep4_D3T12")
K27.rep1.rep2 <- merge(K27.rep1, K27.rep2, 
                      by=c("K27_chr","K27_start","K27_end"), 
                      all=F, sort=F)
K27.rep1.rep2.rep3 <- merge(K27.rep1.rep2, K27.rep3,
                           by=c("K27_chr","K27_start","K27_end"), 
                           all=F, sort=F)
K27.rep1.rep2.rep3.rep4 <- merge(K27.rep1.rep2.rep3, K27.rep4,
                                by=c("K27_chr","K27_start","K27_end"), 
                                all=F, sort=F)

gene <- read.table("Ahal_v2_2_1_gene.bed", sep="\t", header=F)
colnames(gene) <- c("K27_chr","K27_start","K27_end","Ahal_ID",
                    "score","strand","source","feature","frame","group")
gene[,2] <- gene[,2]+1

K27.rep1.rep2.rep3.rep4.gene <- merge(K27.rep1.rep2.rep3.rep4, gene, 
                                     by=c("K27_chr","K27_start","K27_end"), 
                                     all=F, sort=F)
K27.rep1.rep2.rep3.rep4.gene <- 
  K27.rep1.rep2.rep3.rep4.gene[!duplicated(paste(K27.rep1.rep2.rep3.rep4.gene$K27_chr,
                                                K27.rep1.rep2.rep3.rep4.gene$K27_start,
                                                K27.rep1.rep2.rep3.rep4.gene$K27_end,
                                                K27.rep1.rep2.rep3.rep4.gene$Ahal_ID),sep=","),]
K27.rep1.rep2.rep3.rep4.gene <- 
  K27.rep1.rep2.rep3.rep4.gene[,-((ncol(K27.rep1.rep2.rep3.rep4.gene)-5):
                                   ncol(K27.rep1.rep2.rep3.rep4.gene))]

K27.mat <- K27.rep1.rep2.rep3.rep4.gene[,4:35]
row.names(K27.mat) <- K27.rep1.rep2.rep3.rep4.gene[,36]
head(K27.mat)

#group <- factor(rep(c("Nov.6", "Dec.4", "Jan.8", "Feb.5", "Mar.5", "Apr.2", 
#                      "Apr.30", "May.28", "Jul.2", "Jul.30", "Aug.27", "Sep.24"),2))
group <- factor(rep(c("A","B","C","D","E","F","G","H"),4))


# GLM approach
design <- model.matrix(~ group)
d <- DGEList(counts = K27.mat, group = group)

keep <- rowSums(cpm(d)>1) >= 4
d <- d[keep, , keep.lib.sizes=FALSE]
Ahal_ID <- row.names(d)
d.counts <- cbind(Ahal_ID, d$counts)
write.csv(d.counts, file="dielK27_CPMover1.csv", row.names=F)

d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)

lrt <- glmLRT(fit, coef = 2:8)
topTags(lrt)
table <- as.data.frame(topTags(lrt, n = nrow(K27.mat)))
Ahal_ID <- row.names(table)
table2 <- cbind(Ahal_ID, table)
write.csv(table2, file="dielK27_CPMover1_edgeR_result.csv", row.names=F)


# Venn diagram
K4.mat <- read.csv("dielK4_CPMover1_edgeR_result.csv",header=T,sep=",")
K27.mat <- read.csv("dielK27_CPMover1_edgeR_result.csv",header=T,sep=",")

K4_diel_osci <- K4.mat[which(K4.mat$FDR<0.05),]
K27_diel_osci <- K27.mat[which(K27.mat$FDR<0.05),]
K4K27_diel_osci <- merge(K4_diel_osci,K27_diel_osci,by=c("Ahal_ID"),all=F,sort=F)


pdf("../figs/edgeR/venn_K4K27_diel_osci.pdf",height=2.5,width=2.5)
draw.pairwise.venn(
  area1=nrow(K4_diel_osci), area2=nrow(K27_diel_osci), cross.area=nrow(K4K27_diel_osci), 
  cex=1, category=c("K4 diel","K27 diel"), cat.cex=c(1,1), 
  fontfamily="sans", col=rep("transparent",2),
  fill=c(rgb(0.9,0,0),rgb(0,0,0.7)), alpha=rep(0.5,2),
  cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.fontfamily="sans",
  cat.col=c(rgb(0.9,0,0),rgb(0,0,0.7)),scaled=T, margin=0.1
)
dev.off()
