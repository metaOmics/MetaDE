##' Heatmap of selected DE genes 
##' Author: Tianzhou Ma
##' Institution: University of pittsburgh
##' Date: 08/20/2016

##' The \code{heatmap.sig.genes} is a function to draw the Heatmap of DE genes 
##' given a FDR cut point obtained from the Meta-analysis.
##' @title A function to plot the heatmap of DE genes detected at a given FDR 
##' threshold from the Meta-analysis.
##' @param result is the output from MetaDE.
##' @param meta.method is the meta-analysis method used in MetaDE.
##' @param fdr.cut is the FDR cutoff used to select the DE genes.
##' @param color is the color of the heatmap.

##' @return a figure shows the standardized expression levels for the DE genes 
##' detected by meta analysis across studies/datasets.
##' @export
##' @examples
##' \dontrun{
##' meta.method <- 'AW'
##' fdr.cut <- 1e-7
##' pdf('heatmap.test.pdf')
##' heatmap.sig.genes(meta.res.p, meta.method=meta.method,
##'                     fdr.cut=fdr.cut,color="GR")  
##' dev.off()
##' }

heatmap.sig.genes<-function(result,meta.method, fdr.cut,color="GR") {
  
    ci<-match(meta.method,colnames(result$meta.analysis$FDR))
    sig.index<-result$meta.analysis$FDR[,ci]<=fdr.cut   
    sig.index[is.na(sig.index)]<-F
    if (sum(sig.index)==0) {
      stop ("0 signficant genes,there's no genes for plotting")
    }
    cat ("# of genes significant=",sum(sig.index),"\n")
    K<-attr(result$meta.analysis,"nstudy")#number of studies
    N<-attr(result$meta.analysis,"nperstudy") #number of samples in each study
    ni<-attr(result$meta.analysis,"nperlabelperstudy")
    ind.method<-attr(result$meta.analysis,"individual.analysis")
    
    #---get standardized data---# 
    sdat<-label<-NULL
    for (i in 1:K) { 
       tempd<-result$raw.data[[i]][[1]]
       templ<-result$raw.data[[i]][[2]]
        colr<-order(templ)
        templ<-templ[colr]
        tempd<-tempd[,colr]
        label<-c(label,templ)    
        sdat<-cbind(sdat,t(scale(t(tempd)))) # scale each study then combine
    }
    sdat<-t(scale(t(sdat)))#scale across all studies 
    #-----plot genes with fdr < cut----------------#
    #if ("AW"%in%attr(result$meta.analysis,"meta.method")|
    #"AW.OC"%in%attr(result$meta.analysis,"meta.method"))
	if (meta.method=="AW") { 
    forplot<-order.genes.AW(dat=sdat[sig.index,],
                           AW.weight=result$meta.analysis$AW.weight[sig.index,])      
   } else {
     forplot<-order.genes.simple(dat=sdat[sig.index,])
   }
    match.index<-match(row.names(forplot),row.names(sdat))
   #----------plot heatmap of significant genes-------------------#
   attr(forplot,"n")<-N
   attr(forplot,"ni")<-ni
   attr(forplot,"label")<-label
   plot.matrix(forplot,color=color)
   #-----------summarize results for signficant genes only--------#
   stat<-pvalue<-NULL
   print(names(result))
   if (!is.null(result$ind.stat)){ 
    stat<-result$ind.stat[match.index,]
    colnames(stat)<-paste("stat",1:K,sep="")
    pvalue<-result$ind.p[match.index,]
    colnames(pvalue)<-paste("pvalue",1:K,sep="")
   }
    meta.stat<-result$meta.analysis$stat[match.index]
    meta.pvalue<-result$meta.analysis$pval[match.index]
    meta.FDR<-result$meta.analysis$FDR[match.index]
    sig.result<-cbind(stat,pvalue,meta.stat,meta.pvalue,meta.FDR)

    #if ("AW"%in%attr(result$meta.analysis,"meta.method")){
    	if (meta.method=="AW"){
    AW.weight<-result$meta.analysis$AW.weight[match.index,]
    colnames(AW.weight)<-paste("W",1:K,sep="")
    sig.result<-cbind(stat,pvalue,meta.stat,meta.pvalue,meta.FDR,AW.weight)
    }  
   return(sig.result)
}