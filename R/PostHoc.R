##' Post-hoc analysis on AW results
##' Author: Tianzhou Ma
##' Institution: University of pittsburgh
##' Date: 08/29/2016

##' The \code{posthoc.aw} is a function to perform post-hoc analysis on AW 
##' results to determine the overall effect size directionality.
##' @title Post-hoc analysis on AW results.
##' @param result is the output from MetaDE AW method
##' @return a new AW result with additional result (the overall effect size
##' directionality) from the post-hoc analysis.
##' @export
##' @examples
##' \dontrun{
##' posthoc.result <- posthoc.aw(result=meta.res)
##' }


posthoc.aw <- function(result) {
	meta.method <- attr(result$meta.analysis,"meta.method")
   	if(meta.method != "AW") {
   		stop ("The post-hoc analysis only available for AW output")
   	}
   	
   	K<-attr(result$meta.analysis,"nstudy")#number of studies
    N<-attr(result$meta.analysis,"nperstudy") #number of samples in each study
    ni<-attr(result$meta.analysis,"nperlabelperstudy")
    data.type <- attr(result$meta.analysis,"data.type")
    dataset.name <- attr(result$meta.analysis,"dataset.name")
    group.name <- attr(result$meta.analysis,"group.name")
    ref.group <- attr(result$meta.analysis,"ref.group")    
    AW.weight <- result$meta.analysis$AW.weight
    AW.cate <- apply(AW.weight,1,paste,collapse=',')
    ind.stat <- result$ind.stat
    raw.data <- result$raw.data
    G <- nrow(raw.data[[1]][[1]])
    
    out.es <- rep(NA,G)    
    tab.cate <- names(table(AW.cate))
    for (w in tab.cate) {
      w.index <- which(AW.cate==w)
      if(sum(AW.weight[w.index[1],])==1) {
      	  for (i in w.index) {
             wt <- AW.weight[i,]      	
             out.es[i] <- ind.stat[i,which(wt==1)]    
          }
       } else {
          study.index <- which(AW.weight[w.index[1],]==1)
          full_dat<-vector(mode = "list", length = length(study.index))
          N <- n <- c()
       for(i in study.index){
           groupLabel <- raw.data[[i]][[2]]
           y<-raw.data[[i]][[1]][w.index,]
           full_dat[[which(study.index==i)]][[1]] <- y
           full_dat[[which(study.index==i)]][[2]] <- groupLabel
           nns<- get.sample.label.number(groupLabel)
           N<-c(N,nns$N) #sample size per study 
           n<-append(n,nns$n) #sample size per label per study
       }
         
       if(data.type=="continuous") {
         ind.res<-ind.cal.ES(full_dat,paired=rep(FALSE,length(study.index)),
                                   nperm=100)
       } else if (data.type=="discrete") {       
        for(i in 1:length(study.index)){
          temp_dat <- full_dat[[i]][[1]]
          libsize <- colSums(temp_dat)
        for (i in 1:nrow(temp_dat)){
           for (j in 1:ncol(temp_dat)) {
             ## obtain log2 count matrix offset by library size
             temp_dat[i,j] <- log2(temp_dat[i,j]+0.25) - log2(libsize[j])
            }  
          }
          full_dat[[i]][[1]] <- temp_dat
        } 
        ind.res<-ind.cal.ES(full_dat,paired=rep(FALSE,length(study.index)),
                                 nperm=100)
      }
        meta.es<-MetaDE.ES(ind.res,meta.method="REM",REM.type="HO")
        out.es[w.index] <- meta.es$zval
    }       	
  }    
    posthoc.stat <- out.es
    posthoc.dir <- rep(NA, length(posthoc.stat))
    posthoc.dir[sign(posthoc.stat)==1] <- "Upregulated"
    posthoc.dir[sign(posthoc.stat)== -1] <- "Downregulated"
    posthoc.result <- data.frame(posthoc.stat=posthoc.stat,
                                 posthoc.dir =posthoc.dir)
    rownames(posthoc.result) <- rownames(raw.data[[1]][[1]])                       
    result.new <- result
    result.new$meta.analysis$posthoc <- posthoc.result
    #attr(result.new$meta.analysis,"posthoc.es")<- out.es       	
	return(result.new)
}