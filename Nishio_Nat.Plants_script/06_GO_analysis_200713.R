
##### Download GO annotaion #####
# Download the gaf file from GOConsortium,
# extract only Goid and A.thaliana id, and
# save as "tair_200323ver.txt"


##### Addition of GOterm to the reference gene list #####
setwd("data")
GO <- read.table("tair_200323ver.txt",header=T,sep="\t")
GO <- GO[!duplicated(paste(GO$GOid,GO$Atha.gene)),]
write.table(GO, file="tair_200323ver_uniq.txt", sep="\t", quote=F, row.names=F)

gene <- read.csv("K4_rep1_rpkm+1_log2_annotation.csv",header=T,sep=",")
gene2 <- gene[,1:2]
gene.GO <- merge(gene2, GO, by.x="Araport11_ID", by.y="Atha.gene",all=F,sort=F)
write.csv(gene.GO, file="Araport11ID_AhalID_GOid.csv", quote=F, row.names=F)



##### Install the following packages of bioconductor #####
# multtest, GO.db
source("http://bioconductor.org/biocLite.R")
biocLite("multtest", "GO.db")

# load the GO analysis functions
source("../functions/GOanalysis_functions.R")



###### GO analysis of the genes with dual seasonality #####

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with dual seasonality
gl <- 
  read.csv("K4K27_data_stat_maxK4K27over2_ampK4K27over1.csv",
           header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
colnames(gl) <- NULL
gl <- gl[,1]
ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a csv file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="Dualseasonality.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")



###### GO analysis of the genes with dual seasonality with phase difference #####
# relative area over 0.15

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
# Genes with dual seasonality
gl <- 
  read.csv("K4K27_stat_maxK4K27over2_ampK4K27over1_cor_relativearea_dir.csv",
           header=T,sep=",")
ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
gl <- gl[gl$Direction=="Clock",1]

ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.5), 
            file="Dualseasonality.phasedif.Araport11ID.GOid.result.0.5.txt", 
            quote=F, row.names=F, sep="\t")



###### GO analysis of the genes with K27 seasonality #####

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with K27 seasonality
gl <- read.csv("K4K27_data_stat_maxK27over2_ampK27over1.csv",header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
colnames(gl) <- NULL
gl <- gl[,1]

ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="K27seasonality.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")



###### GO analysis of the genes with K27 seasonality (winter and summer peak) #####
#　Winter

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with K27 seasonality
gl <- read.csv("K4K27_monthmean_ampK27over1_peakmonth.csv",header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
gl <- subset(gl, gl$K27.maxK27 == 12 | gl$K27.maxK27 == 1 | gl$K27.maxK27 == 2)
colnames(gl) <- NULL
gl <- gl[,1]

ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="K27seasonality.winter.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")


# Summer

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with K27 seasonality
gl <- read.csv("K4K27_monthmean_ampK27over1_peakmonth.csv",header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
gl <- subset(gl, gl$K27.maxK27 == 6 | gl$K27.maxK27 == 7 | gl$K27.maxK27 == 8)
colnames(gl) <- NULL
gl <- gl[,1]
ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="K27seasonality.summer.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")



###### GO analysis of the genes with K4 seasonality #####
# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with K4 seasonality
gl <- read.csv("K4K27_data_stat_maxK4over2_ampK4over1.csv",header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
colnames(gl) <- NULL
gl <- gl[,1]

ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="K4seasonality.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")



###### GO analysis of the genes with K4 seasonality (winter and summer peak) #####
#　Winter

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with K4 seasonality
gl <- read.csv("K4K27_monthmean_ampK4over1_peakmonth.csv",header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
gl <- subset(gl, gl$K4.maxK4 == 12 | gl$K4.maxK4 == 1 | gl$K4.maxK4 == 2)
colnames(gl) <- NULL
gl <- gl[,1]

ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="K4seasonality.winter.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")


# Summer

# load the table of Gene ID and GO (as ulg)
ulg <- read.csv("Araport11ID_AhalID_GOid.csv",header=T,sep=",")

# load the gene list (as gl)
#  Genes with K4 seasonality
gl <- read.csv("K4K27_monthmean_ampK4over1_peakmonth.csv",header=T,sep=",")

ulg <- ulg[,2:3]
colnames(ulg) <- c("locus","GOid")
gl <- subset(gl, gl$K4.maxK4 == 6 | gl$K4.maxK4 == 7 | gl$K4.maxK4 == 8)
colnames(gl) <- NULL
gl <- gl[,1]
ulg<-as.matrix(ulg)
gl <- as.matrix(gl)

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a txt file
write.table(ng.prepGOtestOutTable(result,0.2), 
            file="K4seasonality.summer.Araport11ID.GOid.result.0.2.txt", 
            quote=F, row.names=F, sep="\t")



##### Visualization of the results #####
dir.create("../figs/GO_barplot")


Dualseasonality.GO <- 
  read.table("Dualseasonality.Araport11ID.GOid.result.0.2.txt",header=T,sep="\t")

pdf("../figs/GO_barplot/Dualseasonality.Araport11ID.GOid.result.0.2.pdf",height=1.2,width=3)
par(mar=c(2.2,10,1.7,1.8))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)
names=c("Response to herbivore",
        "Green leaf volatile biosynthetic process",
        "Response to wounding",
        "Response to abscisic acid"
)
barplot=barplot(-log10(Dualseasonality.GO[c(4:1),1]),xlim=c(0,5),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1,
                space=0.4, horiz=T, border=F, axes=F,col="grey") 
axis(side=1, pos=0, tcl=-0.1, labels=F)
axis(side=1, at=seq(0,5,1), pos=0, tcl=-0.1, labels=F)
mtext(c("0","1","2","3","4","5"),at=seq(0,5,1),side=1,line=-0.25)
axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
mtext("Enrichment score",side=1,line=0.2)
mtext(expression(paste("[-", log[10],"(adjusted ",
                       italic("P"),")]", 
                       sep="")),side=1,line=0.75)
mtext(names,at=barplot,side=2,line=0.2,las=1)
arrows(0, 0, 0, 6, length=0)
mtext("Genes with",side=3,line=0.55)
mtext("Dual seasonality",side=3,line=0.1)
arrows(-log10(0.05),0,-log10(0.05),6, length=0, lty=2, 
       col="black")
dev.off()


Dualseasonality.GO <- 
  read.table("Dualseasonality.phasedif.Araport11ID.GOid.result.0.5.txt",header=T,sep="\t")

pdf("../figs/GO_barplot/Dualseasonality.phasedif.Araport11ID.GOid.result.0.5.pdf",height=1.2,width=3)
par(mar=c(2.2,10,1.7,1.8))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)
names=c("Cold acclimation",
        "Negative regulation of flower development",
        "Green leaf volatile biosynthetic process",
        "Vernalization response"
)
barplot=barplot(-log10(Dualseasonality.GO[c(4:1),1]),xlim=c(0,2.5),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, horiz=T, border=F, axes=F,col="grey") 
axis(side=1, at=c(0,1,2), pos=0, tcl=-0.1, labels=F)
arrows(0, 0, 2.5, 0, length=0)
mtext(c("0","1","2"),at=c(0,1,2),side=1,line=-0.25)
axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
mtext("Enrichment score",side=1,line=0.2)
mtext(expression(paste("[-", log[10],"(adjusted ",
                       italic("P"),")]", 
                       sep="")),side=1,line=0.75)
mtext(names,at=barplot,side=2,line=0.2,las=1)
arrows(0, 0, 0, 6, length=0)
mtext("Genes with",side=3,line=0.55)
mtext("phase differences",side=3,line=0.1)
arrows(-log10(0.05),0,-log10(0.05),6, length=0, lty=2, 
       col="black")
dev.off()


K27seasonality.GO <- 
  read.table("K27seasonality.Araport11ID.GOid.result.0.2.txt",header=T,sep="\t")

pdf("../figs/GO_barplot/K27seasonality.Araport11ID.GOid.result.0.2.pdf",height=0.97,width=3)
par(mar=c(2.2,10,1.7,1.8))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)
names=c("Green leaf volatile biosynthetic process",
        "Response to herbivore"
)
barplot=barplot(-log10(K27seasonality.GO[c(8,4),1]),xlim=c(0,2),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, horiz=T, border=F, axes=F,col=rgb(0,0,0.7)) 
axis(side=1, at=c(0,1,2), pos=0, tcl=-0.1, labels=F)
mtext(c("0","1","2"),at=c(0,1,2),side=1,line=-0.25)
axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
mtext("Enrichment score",side=1,line=0.2)
mtext(expression(paste("[-", log[10],"(adjusted ",
                       italic("P"),")]", 
                       sep="")),side=1,line=0.75)
mtext(names,at=barplot,side=2,line=0.2,las=1)
arrows(0, 0, 0, 3.2, length=0)
mtext("Genes with",side=3,line=0.55)
mtext("H3K27me3 seasonality",side=3,line=0.1)
arrows(-log10(0.05),0,-log10(0.05),3.2, length=0, lty=2, 
       col="black")
dev.off()


K27seasonality.GO <- 
  read.table("K27seasonality.winter.Araport11ID.GOid.result.0.2.txt",header=T,sep="\t")

pdf("../figs/GO_barplot/K27seasonality.winter.Araport11ID.GOid.result.0.2.pdf",height=1.2,width=3)
par(mar=c(2.2,10,1.7,1.8))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)
names=c("Response to jasmonic acid",
        "Response to wounding",
        "Response to herbivore",
        "Green leaf volatile biosynthetic process"
)
barplot=barplot(-log10(K27seasonality.GO[c(6,4:2),1]),xlim=c(0,2.5),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, horiz=T, border=F, axes=F,col=rgb(0,0,0.7)) 
axis(side=1, at=c(0,1,2), pos=0, tcl=-0.1, labels=F)
mtext(c("0","1","2"),at=c(0,1,2),side=1,line=-0.25)
arrows(0, 0, 2.5, 0, length=0)
axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
mtext("Enrichment score",side=1,line=0.2)
mtext(expression(paste("[-", log[10],"(adjusted ",
                       italic("P"),")]", 
                       sep="")),side=1,line=0.75)
mtext(names,at=barplot,side=2,line=0.2,las=1)
arrows(0, 0, 0, 6, length=0)
mtext("Genes with",side=3,line=0.55)
mtext("H3K27me3 seasonality (winter)",side=3,line=0.1)
arrows(-log10(0.05),0,-log10(0.05),6, length=0, lty=2, 
       col="black")
dev.off()


K4seasonality.GO <- 
  read.table("K4seasonality.winter.Araport11ID.GOid.result.0.2.txt",header=T,sep="\t")

pdf("../figs/GO_barplot/K4seasonality.winter.Araport11ID.GOid.result.0.2.pdf",height=1.55,width=3)
par(mar=c(2.2,10,1.7,1.8))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)
names=c("DNA unwinding involved in DNA replication",
        "Response to flooding",
        "DNA replication initiation",
        "Resolution of meiotic recombination intermediates",
        "Double-strand break repair via homologous recombination",
        "Reciprocal meiotic recombination",
        "Response to abscisic acid"
)
barplot=barplot(-log10(K4seasonality.GO[c(13,12,11,9,7,5,1),1]),xlim=c(0,6.5),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, horiz=T, border=F, axes=F,col=rgb(0.9,0,0)) 
axis(side=1, at=seq(0,6,2), pos=0, tcl=-0.1, labels=F)
mtext(c("0","2","4","6"),at=seq(0,6,2),side=1,line=-0.3)
arrows(0, 0, 6.5, 0, length=0)
axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
mtext("Enrichment score",side=1,line=0.15)
mtext(expression(paste("[-", log[10],"(adjusted ",
                       italic("P"),")]", 
                       sep="")),side=1,line=0.7)
mtext(names,at=barplot,side=2,line=0.2,las=1)
arrows(0, 0, 0, 10.2, length=0)
mtext("Genes with",side=3,line=0.55)
mtext("H3K4me3 seasonality (winter)",side=3,line=0.1)
arrows(-log10(0.05),0,-log10(0.05),10.2, length=0, lty=2, 
       col="black")
dev.off()


K4seasonality.GO <- 
  read.table("K4seasonality.summer.Araport11ID.GOid.result.0.2.txt",header=T,sep="\t")

pdf("../figs/GO_barplot/K4seasonality.summer.Araport11ID.GOid.result.0.2.pdf",height=3,width=3)
par(mar=c(2.2,10,1.7,1.8))
par(mgp=c(0.7,0.2,0.1))
par(ps=6)
par(cex=1)
par(xpd=T)
names=c("Response to oomycetes",
        "Regulation of gibberellin biosynthetic process",
        "Salicylic acid mediated signaling pathway",
        "Response to salicylic acid",
        "Response to water deprivation",
        "Response to bacterium",
        "Regulation of jasmonic acid mediated signaling pathway",
        "Regulation of defense response",
        "Response to chitin",
        "Response to heat",
        "Defense response to fungus",
        "Response to reactive oxygen species",
        "Response to salt stress",
        "Regulation of systemic acquired resistance",
        "Response to herbivore",
        "Protein complex oligomerization",
        "Response to hydrogen peroxide",
        "Response to jasmonic acid",
        "Cellular response to hypoxia",
        "Response to wounding"
)
barplot=barplot(-log10(K4seasonality.GO[c(24,23,21,20,18,16,14:1),1]),xlim=c(0,9.5),
                names.arg=F, xlab=NULL, las=1, tcl=-0.1, 
                space=0.4, horiz=T, border=F, axes=F,col=rgb(0.9,0,0)) 
axis(side=1, at=seq(0,9,3), pos=0, tcl=-0.1, labels=F)
mtext(c("0","3","6","9"),at=seq(0,9,3),side=1,line=-0.57)
arrows(0, 0, 9.5, 0, length=0)
axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
mtext("Enrichment score",side=1,line=-0.12)
mtext(expression(paste("[-", log[10],"(adjusted ",
                       italic("P"),")]", 
                       sep="")),side=1,line=0.43)
mtext(names,at=barplot,side=2,line=0.2,las=1)
arrows(0, 0, 0, 28.5, length=0)
mtext("Genes with",side=3,line=0.33)
mtext("H3K4me3 seasonality (summer)",side=3,line=-0.12)
arrows(-log10(0.05),0,-log10(0.05),28.5, length=0, lty=2, 
       col="black")
dev.off()
