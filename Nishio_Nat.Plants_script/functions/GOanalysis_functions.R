
################################
##### GOanalysis_functions #####
################################

library(multtest)
library(GO.db)


# ng.BHFDR ------------
# A function to adjust p-value for multiple comparison
# Benjamini & Hochberg (1995) step-up FDR
# Out put adjusted p-value after correction of FDR 

ng.BHFDR <- function(rawp){
  tmp <- mt.rawp2adjp(rawp, proc=c("BH"))
  adjp <- tmp$adjp
  index <- tmp$index
  out <- adjp[order(index), 2]
  return(out)
}



# multiple fisher test --------------------------
ng.mft <- function(
  cgt, #output from ng.MakeContigGOidTable, [,1]:"locus", [,2]:"GOid"
  gn.test, #contig names for test
  alternative="greater"
){
  
  #cat(sprintf("%s\n", Sys.time()))
  
  gid.u <- unique(cgt[,"GOid"])
 
  ft.in <- matrix(0, nrow=length(gid.u), ncol=9)
  colnames(ft.in) <- c("xtt", "xft", "xtf", "xff", "xnt", "xnf", "xtn", "xfn", "xnn")
  rownames(ft.in) <- gid.u
  
  #           　　gn.test
  #       　　　 TRUE FALSE
  #Group　TRUE   xtt   xft   xnt
  #       FALSE  xtf   xff   xnf
  #              xtn   xfn   xnn

  ft.in[,"xnn"] <- length(unique(cgt[, "locus"]))
  gn.pp.gid <- table(cgt[, "GOid"])
  ft.in[names(gn.pp.gid), "xnt"] <- gn.pp.gid
  ft.in[,"xnf"] <- ft.in[,"xnn"] - ft.in[,"xnt"]
  ft.in[,"xtn"] <- length(intersect(gn.test, unique(cgt[, "locus"])))
  ft.in[,"xfn"] <- ft.in[,"xnn"] - ft.in[,"xtn"]
  
  gsea.test <- cgt[is.element(cgt[,"locus"], gn.test), ]
  gn.test.gid <- table(gsea.test[, "GOid"])
  ft.in[names(gn.test.gid), "xtt"] <- gn.test.gid

  ft.in[,"xtf"] <- ft.in[,"xtn"] - ft.in[,"xtt"]
  ft.in[,"xft"] <- ft.in[,"xnt"] - ft.in[,"xtt"]
  ft.in[,"xff"] <- ft.in[,"xnf"] - ft.in[,"xtf"]
  
  #Fisher's exact test.  8? sec
  fr <- rep(1, nrow(ft.in))
  dt <- rep(1, nrow(ft.in))
  for(i in 1:nrow(ft.in)){
    start <- Sys.time()
    if(ft.in[i,"xtn"] > 1 && ft.in[i,"xnt"] > 1){ 
      contable <- matrix(ft.in[i, 1:4], ncol=2)
      tmp <- fisher.test(contable, alternative = alternative)
      fr[i] <- tmp$p.value
    } else {
    }
    end <- Sys.time()
    dt[i] <- end - start
  }
  
  out <- cbind(fr, ft.in, dt)
  colnames(out) <- c("p.value", colnames(ft.in), "time")
  rownames(out) <- rownames(ft.in)
  
  return(out)
  
}


# get GO terms -----------------------------
ng.GetGOTerms <- function(GOid){
  
  out <- matrix(NA, nrow=length(GOid), ncol=4)
  for(i in 1:length(GOid)){
    tmp <- try(get(GOid[i], GOTERM), silent=TRUE)
    if(class(tmp)=="try-error"){
      out[i,] <- rep("NA",length=4)
    } else {
      out[i,] <- c(GOid[i], Term(tmp), Ontology(tmp), Definition(tmp))
    }
  }
  return(out)
}

#out <- c(out, tmp, Term(tmp), Synonym(tmp), Definition(tmp), Ontology(tmp))



# prep. GO test output table ------------------------------

ng.prepGOtestOutTable <- function(r, alpha){
  adp <- ng.BHFDR(r[,"p.value"])
  tmp.id <- rownames(r)[adp<alpha]
  tmp.adp <- adp[adp<alpha]
  if(length(tmp.id)!=0){
  tmp.description <- ng.GetGOTerms(tmp.id)
  }
  tmp.xnn <- r[adp<alpha, c("xtt", "xtn", "xnt", "xnn")]
  if(length(tmp.id)==0){
    out <- "No GO term enriched"
  }else if(length(tmp.id)==1){
    out <- c(tmp.adp, tmp.id, c(tmp.description), tmp.xnn)
    out <- as.matrix(out)
    row.names(out) <- c("Adjusted_P_value", "ID", "ID2", "Term", "Ontology", 
                        "Definition", "A & B", "A", "B", "U")
    out <- t(out)
    out <- na.omit(out)
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
    out <- as.data.frame(out)
    sortlist <- order(out$tmp.adp)
    out <- out[sortlist,]
    colnames(out) <- c("Adjusted_P_value", "ID", "ID2", "Term", "Ontology", 
                       "Definition", "A & B", "A", "B", "U")
    out <- na.omit(out)
  }
  return(out)
}



# Function to transform the first character to uppercase
first.upper <- function(cha){
	h <- substr(cha,1,1)
	t <- substr(cha,2,nchar(cha))
	out <- paste0(toupper(h),t)
	return(out)
}



# Plot the result
plot.func <- function(filename_wo_extension, mod){

  GO <- read.table(paste0(filename_wo_extension,".txt"),header=T,sep="\t",quote="")

  if(mod=="H3K4me3"){
    col<-rgb(0.55,0,0)
  }else if(mod=="H3K27me3"){
    col<-rgb(0,0,1)
  }else{
    col<-rgb(0,0,0.5)
  }
  
if(length(GO$Term[GO$Ontology=="BP"])>=4){
  pdf(paste0("../figs/GO_barplot/",filename_wo_extension,".pdf"),height=1.2,width=3)
  par(mar=c(2.2,10,1.7,1.8))
  par(mgp=c(0.7,-0.32,0.1))
  par(ps=6)
  par(cex=1)
  par(xpd=T)

  names <- GO$Term[GO$Ontology=="BP"][4:1]
  max.val <- ceiling(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][1]))

  barplot=barplot(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][4:1]),
                  xlim=c(0,max.val),
                  names.arg=F, xlab=NULL, las=1, tcl=-0.1,
                  space=0.4, horiz=T, border=F, axes=F,col=col) 
  axis(side=1, pos=0, tcl=-0.1)
  #mtext(as.character(0:max.val),at=seq(0,max.val,1),side=1,line=-0.25)
  axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
  mtext("Enrichment score",side=1,line=0.2)
  mtext(expression(paste("(-", log[10],"[Adjusted ",
                         italic("P")," value])", 
                         sep="")),side=1,line=0.75)
  for(j in 1:length(names)){
  	mtext(first.upper(as.character(names[j])),at=barplot[j],side=2,line=0.2,las=1)
  }
  arrows(0, 0, 0, 6, length=0)
  mtext(paste0("Category ",i),side=3,line=0.1)
  arrows(-log10(0.05),0,-log10(0.05),6, length=0, lty=2, 
       col="black")
  dev.off()
}else if(length(GO$Term[GO$Ontology=="BP"])==3){
  pdf(paste0("../figs/GO_barplot/",filename_wo_extension,".pdf"),height=1.1,width=3)
  par(mar=c(2.2,10,1.7,1.8))
  par(mgp=c(0.7,-0.32,0.1))
  par(ps=6)
  par(cex=1)
  par(xpd=T)
  
  names <- GO$Term[GO$Ontology=="BP"][3:1]
  max.val <- ceiling(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][1]))
  
  barplot=barplot(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][3:1]),
                  xlim=c(0,max.val),
                  names.arg=F, xlab=NULL, las=1, tcl=-0.1,
                  space=0.4, horiz=T, border=F, axes=F,col=col) 
  axis(side=1, pos=0, tcl=-0.1)
  #mtext(as.character(0:max.val),at=seq(0,max.val,1),side=1,line=-0.25)
  axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
  mtext("Enrichment score",side=1,line=0.2)
  mtext(expression(paste("(-", log[10],"[Adjusted ",
                         italic("P")," value])", 
                         sep="")),side=1,line=0.75)
  for(j in 1:length(names)){
  	mtext(first.upper(as.character(names[j])),at=barplot[j],side=2,line=0.2,las=1)
  }
  arrows(0, 0, 0, 4.5, length=0)
  mtext(paste0("Category ",i),side=3,line=0.1)
  arrows(-log10(0.05),0,-log10(0.05),4.5, length=0, lty=2, 
         col="black")
  dev.off()
}else if(length(GO$Term[GO$Ontology=="BP"])==2){
  pdf(paste0("../figs/GO_barplot/",filename_wo_extension,".pdf"),height=0.975,width=3)
  par(mar=c(2.2,10,1.7,1.8))
  par(mgp=c(0.7,-0.32,0.1))
  par(ps=6)
  par(cex=1)
  par(xpd=T)
  
  names <- GO$Term[GO$Ontology=="BP"][2:1]
  max.val <- ceiling(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][1]))
  
  barplot=barplot(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][2:1]),
                  xlim=c(0,max.val),
                  names.arg=F, xlab=NULL, las=1, tcl=-0.1,
                  space=0.4, horiz=T, border=F, axes=F,col=col) 
  axis(side=1, pos=0, tcl=-0.1)
  #mtext(as.character(0:max.val),at=seq(0,max.val,1),side=1,line=-0.25)
  axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
  mtext("Enrichment score",side=1,line=0.2)
  mtext(expression(paste("(-", log[10],"[Adjusted ",
                         italic("P")," value])", 
                         sep="")),side=1,line=0.75)
  for(j in 1:length(names)){
  	mtext(first.upper(as.character(names[j])),at=barplot[j],side=2,line=0.2,las=1)
  }
  arrows(0, 0, 0, 3.1, length=0)
  mtext(paste0("Category ",i),side=3,line=0.1)
  arrows(-log10(0.05),0,-log10(0.05),3.1, length=0, lty=2, 
         col="black")
  dev.off()
}else if(length(GO$Term[GO$Ontology=="BP"])==1){
  pdf(paste0("../figs/GO_barplot/",filename_wo_extension,".pdf"),height=0.85,width=3)
  par(mar=c(2.2,10,1.7,1.8))
  par(mgp=c(0.7,-0.32,0.1))
  par(ps=6)
  par(cex=1)
  par(xpd=T)
  
  names <- GO$Term[GO$Ontology=="BP"][1]
  max.val <- ceiling(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][1]))
  
  barplot=barplot(-log10(GO$Adjusted_P_value[GO$Ontology=="BP"][1]),
                  xlim=c(0,max.val),
                  names.arg=F, xlab=NULL, las=1, tcl=-0.1,
                  space=0.4, horiz=T, border=F, axes=F,col=col) 
  axis(side=1, pos=0, tcl=-0.1)
  #mtext(as.character(0:max.val),at=seq(0,max.val,1),side=1,line=-0.25)
  axis(side=2, pos=0, tcl=-0.1, label=F, at=barplot)
  mtext("Enrichment score",side=1,line=0.2)
  mtext(expression(paste("(-", log[10],"[Adjusted ",
                         italic("P")," value])", 
                         sep="")),side=1,line=0.75)
  for(j in 1:length(names)){
  	mtext(first.upper(as.character(names[j])),at=barplot[j],side=2,line=0.2,las=1)
  }
  arrows(0, 0, 0, 1.5, length=0)
  mtext(paste0("Category ",i),side=3,line=0.1)
  arrows(-log10(0.05),0,-log10(0.05),1.5, length=0, lty=2, 
         col="black")
  dev.off()
}
  
}