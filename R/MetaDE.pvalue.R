##' Meta analysis by combining p-value 
##' The \code{MetaDE} is a function to identify genes associated with the 
##' response/phenoype of interest (can be either group, continuous or survival) 
##' by combining p-values from multiple studies(datasets). 
##' The main input consists of p-values from your own method/calculations.
##' @title Meta analysis by combining p-value
##' @param x is a list with components:
##' \itemize{
##' \item{p:}{ a list of p values for each dataset.}
##' \item{bp:}{  a list of p values calculated from permutation for each dataset. 
##' This part can be NULL if you just have the p-values from your own method.}
##' }
##' @param x is a list with components:
##' @param meta.method is a character to specify the Meta-analysis method 
##' used to combine the p-values.
##' @param rth is the option for roP and roP.OC method. rth means the 
##' rth smallest p-value. 
##' @param parametric is a logical values indicating whether the parametric 
##' methods is chosen to calculate the p-values in meta-analysis.
##' @return a list with components: \cr
##' \itemize{
##' \item{stat:}{ a matrix with rows representing genes. It is the statistic for 
##' the selected meta analysis method of combining p-values.}
##' \item{pval:}{ the p-value from meta analysis for each gene for the above 
##' stat.}
##' \item{FDR:}{ the FDR of the p-value for each gene for the above stat.}
##' \item{AW.weight:}{ The optimal weight assigned to each dataset/study for 
##' each gene if the '\code{AW}' method was chosen.}
##' }
##' @export
##' @examples
##' data('Leukemia')
##' data('LeukemiaLabel')
##' data <- Leukemia
##' K <- length(data)
##' clin.data <- lapply(label, function(x) {data.frame(x)} )
##' for (k in 1:length(clin.data)){
##'  colnames(clin.data[[k]]) <- "label"
##' }
##' select.group <- c('inv(16)','t(15;17)')
##' ref.level <- "inv(16)"
##' data.type <- "continuous"
##' ind.method <- c('limma','limma','limma')
##' resp.type <- "twoclass"
##' paired <- rep(FALSE,length(data))
##' ind.res <- Indi.DE.Analysis(data=data,clin.data= clin.data, 
##'                         data.type=data.type,resp.type = resp.type,
##'                         response='label',
##'                         ind.method=ind.method,select.group = select.group,
##'                         ref.level=ref.level,paired=paired)
##' meta.method <- "AW"
##' meta.res <- MetaDE.pvalue(ind.res,meta.method,rth=NULL,parametric=TRUE)
##' summary <- data.frame(ind.p = meta.res$ind.p,
##'                       stat = meta.res$meta.analysis$stat,
##'                       pval = meta.res$meta.analysis$pval,
##'                       FDR = meta.res$meta.analysis$FDR,
##'                       weight = meta.res$meta.analysis$AW.weight)
                       
MetaDE.pvalue <-function(x,meta.method,rth=NULL,parametric=TRUE) {
  #meta.method<-match.arg(meta.method,several.ok = TRUE)
  check.parametric(meta.method,parametric)
  K<-ncol(x$p)
  if (parametric) x$bp<-NULL     
  nm<-length(meta.method)
  meta.res<-list(stat=NA,pval=NA,FDR=NA,AW.weight=NA)
  meta.res$stat<-meta.res$pval<-meta.res$FDR<-matrix(NA,nrow(x$p),nm)
  for( i in 1:nm){
    temp<-switch(meta.method[i],
                 maxP={get.maxP(x$p,x$bp)},minP={get.minP(x$p,x$bp)},
                 Fisher={get.fisher(x$p,x$bp)},roP={get.roP(x$p,x$bp,rth=rth)},
                 AW={get.AW(x$p)},
                 Fisher.OC={get.fisher.OC(x$p,x$bp)},
                 maxP.OC={get.maxP.OC(x$p,x$bp)},
                 minP.OC={get.minP.OC(x$p,x$bp)},
                 roP.OC={get.roP.OC(x$p,x$bp,rth=rth)},
                 Stouffer={get.Stouff(x$p,x$bp)},
                 Stouffer.OC={get.Stouff.OC(x$p,x$bp)},
                 SR={get.SR(x$p,x$bp)},PR={get.PR(x$p,x$bp)})
    meta.res$stat[,i]<-temp$stat
    meta.res$pval[,i]<-temp$pval
    meta.res$FDR[,i]<-temp$FDR
    if(meta.method[i]=="AW"){
      meta.res$AW.weight<-temp$AW.weight
    }
  }
  colnames(meta.res$stat)<-colnames(meta.res$pval)<-colnames(meta.res$FDR)<-
    meta.method
  rownames(meta.res$stat)<-rownames(meta.res$pval)<-
    rownames(meta.res$FDR)<-rownames(x$p)   
  attr(meta.res,"nstudy")<-K
  attr(meta.res,"meta.method")<-meta.method 
  res<-list(meta.analysis=meta.res,ind.p=x$p)	 
  #class(res)<-"MetaDE.pvalue"
  return(res)
}
