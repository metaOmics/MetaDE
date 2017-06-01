##' Function to summarize the results in tabular form 
##' The \code{summary.meta} is a function to summarize the meta-analysis 
##' results from the MetaDE output.
##' @title Function to summarize the meta-analysis results.
##' @param result is the output from MetaDE
##' @param meta.method is the meta-analysis method used in MetaDE
##' @param resp.type is a character indicating the response type, 
##' must be one of "twoclass", "multiclass", "continuous" and "survival".
##' @return a summary table including individual study statistics and pvalue,
##' meta-analysis test statistics, pvalue, FDR, AW weights etc.  
##' @export
##' @examples
##' \dontrun{
##' meta.method <- 'AW'
##' resp.type <- "twoclass"
##' summary.result <- summary.meta(result=meta.res, 
##'                               meta.method = meta.method,
##'                               resp.type = resp.type)
##' posthoc.result <- posthoc.aw(result=meta.res)
##' summary.posthoc.result <- summary.meta(result=posthoc.result, 
##'                               meta.method = meta.method,
##'                               resp.type = resp.type)
##' }

summary.meta <- function(result,meta.method,resp.type) {
	
  if (resp.type!="twoclass" && meta.method != "minMCC"){
  	    summary <- data.frame(ind.p = result$ind.p,
                          pval = result$meta.analysis$pval,
                          FDR = result$meta.analysis$FDR)
    colnames(summary)[(ncol(summary)-1):(ncol(summary))] <- 
      c("pval","FDR")
  }	else {
  	
  if ("minMCC"%in%meta.method) {
    summary <- data.frame(stat = result$meta.analysis$stat,
                          pval = result$meta.analysis$pval,
                          FDR = result$meta.analysis$FDR)
  }
  
  if ("rankProd"%in%meta.method) {
    summary <- data.frame(meta.stat.up = result$meta.analysis$meta.stat.up,
                          pval.up = result$meta.analysis$pval.up,
                          FDR.up = result$meta.analysis$FDR.up,
                          meta.stat.down = result$meta.analysis$meta.stat.down,
                          pval.down = result$meta.analysis$pval.down,
                          FDR.down = result$meta.analysis$FDR.down,
                          AveFC = result$meta.analysis$AveFC)
  }
  
  if ("FEM"%in%meta.method|"REM"%in%meta.method) {
    summary <- data.frame(ind.ES = result$ind.ES, ind.Var=result$ind.Var,
                          zval = result$meta.analysis$zval,
                          pval = result$meta.analysis$pval,
                          FDR = result$meta.analysis$FDR)
    colnames(summary)[(ncol(summary)-2):(ncol(summary))] <- 
      c("zval","pval","FDR")
  }
  
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))==0)  {
    
    summary <- data.frame(ind.stat = result$ind.stat,
                          ind.p = result$ind.p,
                          stat = result$meta.analysis$stat,
                          pval = result$meta.analysis$pval,
                          FDR = result$meta.analysis$FDR)
    colnames(summary)[(ncol(summary)-2):(ncol(summary))] <- 
      c("stat","pval","FDR")
    if (meta.method == "AW") {
      weight <- result$meta.analysis$AW.weight
      colnames(weight) <- paste("weight.study",1:length(result$raw.data),sep="")
      summary <- cbind(summary,weight)
      if(!is.null(result$meta.analysis$posthoc)) {
        posthoc.stat <- result$meta.analysis$posthoc[,1]
        posthoc.dir <- result$meta.analysis$posthoc[,2]
      	summary <- cbind(summary,posthoc.stat,posthoc.dir)
       }
    } 
  }
  
 }
 
  rownames(summary) <- rownames(result$raw.data[[1]][[1]]) 
  return(summary)
}

