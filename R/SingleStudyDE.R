##' Main Function for Individual Study DE: microarray & RNAseq 
##' The \code{Indi.DE.Analysis} is a function to perform individual association
##' analysis between gene expression and the response/phenoype of interest (can 
##' be either group, continuous or survival).
##' @title Main Function for Individual Study DE: microarray & RNAseq. 
##' @param data is a list of K elements, where K is the number of studies, each 
##' element is a microarray or RNAseq expression matrix with G rows and N 
##' columns, where G is number of matched genes and N is the sample size.
##' @param clin.data is a list of K elements, each element includes is a
##' clinical data frame with N rows and p columns, where N is the sample size
##' and p is the number of clinical variables (main response included). 
##' @param data.type is a character indicating the data type of the elements 
##' in \code{data}, must be "continuous" or "discrete".
##' @param resp.type is a character indicating the response type of the 
##' \code{response} variable selected, must be one of "twoclass", "multiclass", 
##' "continuous" and "survival".
##' @param response is one column name of \code{clin.data}, indicating the 
##' phenotype of interest. For survival, two column names have to be specified,
##' the first is the survival time and the second is the censoring status.
##' @param covariate are the clinical covariates you wish to adjust for in the
##' DE analysis, can be a vector of column names or NULL.
##' @param ind.method is a character vector to specify the method used to test 
##' if there is association between the gene expression and outcome variable.
##' must be one of "limma", "sam" for "continuous" data type and "edgeR", 
##' "DESeq2" or "limmaVoom" for "discrete" data type.
##' @param select.group: for two-class comparison only, specify the two groups 
##' for comparison when the group factor has more than two levels.
##' @param ref.level: for two-class/multi-class comparison only, specify the 
##' reference level of the group factor. 
##' @param paired: logical value indicating whether paired design;
##' @param asymptotic: a logical value indicating whether asymptotic 
##' distribution should be used. If FALSE, permutation will be performed. 
## 'For resp.type = "continuous", "survival" only.
##' @param nperm: the number of permutations. Applicable when \code{asymptotic} 
##' is FALSE.
##' @param tail: a character string specifying the alternative hypothesis, 
##' must be one of "abs" (default), "low" or "high". For resp.type = 
##' "continuous", "survival" only.
##' @param seed: Optional initial seed for random number generator.  

##' @return a list with components: \cr
##' \itemize{
##' \item{p}{ For all types of response, the p-value of the association test for 
##' each gene}
##' \item{stat}{ For "continuous" and "survival" only, the value of test 
##' statistic for each ##' gene}
##' \item{bp}{ For "continuous" and "survival" only, the p-value from 
##'	\code{nperm} permutations for each gene. It will be used for the meta 
##' analysis by default. It can be NULL if you chose asymptotic results. }
##' \item{log2FC}{ For "twoclass" only, the log2 fold change for each gene}
##' \item{lfcSE}{ For "twoclass" only, the standard error of log2 fold change 
##' for each gene}
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
##' N <- sapply(data, FUN=function(x) ncol(x))
##' survival.time <- lapply(N,FUN = function(x) round(runif(x,10,2000)))
##' censor.status <- lapply(N,FUN = function(x) sample(c(0,1),x,replace=TRUE) )
##' for (k in 1:length(clin.data)){
##'   clin.data[[k]] <- cbind(clin.data[[k]],survival.time[[k]],censor.status[[k]])
##'   colnames(clin.data[[k]])[2:3] <- c("survival","censor")
##' }
##' ind.method <- c('logrank','logrank','logrank')
##' resp.type <- "survival"
##' ind.res <- Indi.DE.Analysis(data=data,clin.data= clin.data, 
##'                          data.type=data.type,resp.type = resp.type,
##'                          response=c("survival","censor"),
##'                          ind.method=ind.method,asymptotic=TRUE)


Indi.DE.Analysis <- function(data, clin.data, data.type, resp.type, 
                             response,covariate=NULL,ind.method, 
                             select.group=NULL, ref.level=NULL, 
                             paired=NULL,asymptotic=NULL,nperm=NULL, 
                             tail="abs", seed=12345,
                             ... ) {

## Call the packages required 

#library(survival)
#library(limma)
#library(samr)
#library(edgeR)
#library(DESeq2)
               
set.seed(seed)
              	
##Input  
##data: list of data matrices, genes on the rows, samples on the columns;
##clin.data: list of clinical data frames, samples on the rows, clinical
##covariates on the columns;
##data.type = c("continuous","discrete");
##resp.type = c("twoclass", "multiclass", "continuous", "survival");
##response: a character (column name of clin.data) indicating the phenotype 
##of interest (survival (2 column names));
##ind.method = c("limma","sam","limmaVoom","edgeR","DESeq2","pearsonr",
##"spearmanr","logrank");
##select.group: specify the two group names for comparison;
##ref.level: specify the reference level of the group factor;
##paired: logical indicating whether paired design;
##asymptotic: logical whether asymptotic dist should be used;
##nperm: number of permutations;
##tail = c("low", "high", "abs");
  
    
  K<-length(data) # number of studies
  G<-nrow(data[[1]]) # number of matched genes
  response.list <- covariate.list <- vector('list',K)
  if(!is.null(covariate)) covariate.list <-NULL
  for (k in 1:K) {
  	response.list[[k]] <- clin.data[[k]][,response]
  	if(!is.null(covariate)) {
  	 covariate.list[[k]] <-	clin.data[[k]][,covariate]
    }
  }  
##check the individual method for the corresponding resp.type and data.type,
##in addition, check the consistency between resp.type and response.

  check.indmethod(response.list, resp.type, data.type, ind.method, 
                  select.group, tail, paired) 

##verify the correct specification of individual method
  ind.method<- match.arg(ind.method,c("limma","sam","limmaVoom","edgeR",
                         "DESeq2","pearsonr","spearmanr","logrank"),
                         several.ok=TRUE)

  data <-check.exp(data)  # check the gene names for expression data

##Association with groups (DE, ANOVA, etc.)  
 if(resp.type %in% c("twoclass","multiclass") )  {
  log2FC<- lfcSE <-  pvalue <- NULL #prespecify the output matrix
  N<-n<-NULL
   for (k in 1:K) { 
  	 groupLabel <- as.factor(response.list[[k]])
  	 groupName=levels(groupLabel)	
  	 name <- select.group
  	 l<- groupLabel[groupLabel %in% name]
     l<- droplevels(l)
     l <- relevel(l,ref=ref.level)
     y<-data[[k]][,groupLabel %in% name]
     if(is.null(covariate.list[[k]])) {
         c<- NULL
       }      
     if(!is.null(covariate.list[[k]]) && is.vector(covariate.list[[k]]))  {
          c <- covariate.list[[k]][groupLabel %in% name]
       } 
     if(!is.null(covariate.list[[k]]) && !is.vector(covariate.list[[k]])) {
          c <- covariate.list[[k]][groupLabel %in% name,]
      } 
    
 ANOVA<- (resp.type == "multiclass") #indicate whether ANOVA should be performed          
 ind.res<-switch(ind.method[k],
                limma={get.limma(y=y,l=l,c=c,name=name, 
                	               ANOVA=ANOVA,tail=tail,paired=paired[k])},
                sam={get.sam(y=y,l=l,name=name,seed=seed,
                	            ANOVA=ANOVA,tail=tail,paired=paired[k])},
                limmaVoom={get.limmaVoom(y=y,l=l,c=c,name=name,
                	         ANOVA=ANOVA,tail=tail,paired=paired[k])},
                edgeR={get.edgeR(y=y,l=l,c=c,name=name,
                	       ANOVA=ANOVA,tail=tail,paired=paired[k])},
                DESeq2={get.DESeq2(y=y,l=l,c=c,name=name,
                	    ANOVA=ANOVA,tail=tail,paired=paired[k])})

  pvalue<- cbind(pvalue,ind.res$pvalue) #p value

 if(!ANOVA) {
   if(any(grepl('log2FC',names(ind.res))) ) {
     log2FC<- cbind(log2FC,ind.res$log2FC) #log2 fold change 
   } else {
     log2FC<- cbind(log2FC,rep(NA,nrow(y)))
   }  

   if(any(grepl('lfcSE',names(ind.res))) )  {
     lfcSE <- cbind(lfcSE,ind.res$lfcSE) #standard error of log2FC 
    } else {
      lfcSE <- cbind(lfcSE,rep(NA,nrow(y)))
    }  
  }

   cat("dataset",k,"is done\n")
   nns<-get.sample.label.number(l)
	 N<-c(N,nns$N) #sample size per study 
	 n<-append(n,nns$n) #sample size per label per study
} # end of K study for loop
  
  if(is.null(names(data))) {
    colnames(log2FC) <-colnames(lfcSE)<-colnames(pvalue)<-
      paste("dataset",1:K,sep="") # naming studies
	} else {
    colnames(log2FC) <-colnames(lfcSE)<-colnames(pvalue)<-names(data) 
	}

  if(is.null(rownames(data[[1]])) ) {
    rownames(log2FC) <- rownames(lfcSE)<-rownames(pvalue)<-
      paste("gene",1:G,sep="") # naming genes
  } else {
    rownames(log2FC) <-rownames(lfcSE)<-rownames(pvalue)<-rownames(data[[1]]) 
  }

  if(!ANOVA) {
    all.res<-list(log2FC=log2FC,lfcSE=lfcSE,p=pvalue) 
  } else {
    all.res<-list(p=pvalue)
  }
  attr(all.res$p,"nperstudy")<-N
  attr(all.res$p,"nperlabelperstudy")<-n 
  attr(all.res$p,"data.type")<- data.type
  attr(all.res$p,"individual.analysis")<-ind.method
  attr(all.res$p,"response.type")<- resp.type
  attr(all.res$p,"alter.hypothesis")<- tail
  return(all.res)
} #end of twoclass, multiclass comparison   
  
 if(resp.type == c("continuous") )  {
  P<-stat<-P.perm<-NULL #prespecify the output matrix
  N<- NULL
  for (k in 1:K) {   
   y<- data[[k]] 
   r<- response.list[[k]]  
   ind.res<-switch(ind.method[k],
                   pearsonr={get.p.pearsonr(y,r,asymptotic=asymptotic, 
                   	         nperm=nperm,tail=tail)},
                   spearmanr={get.p.spearmanr(y,r,asymptotic= asymptotic,
                   	         nperm=nperm,tail=tail)})  
   P<- cbind(P,(ind.res$obs)[,2]) #p value
   stat<-cbind(stat,(ind.res$obs[,1])) # observed stats
   if (!is.null(nperm)) {
     P.perm<-cbind(P.perm,c(ind.res$perm)) #p value from permutation 
   }
   cat("dataset",k,"is done\n")
	 N<-c(N,ncol(y))	  
} # end of K study for loop
   
  if(is.null(names(data))) {
    colnames(stat)<-colnames(P)<-paste("dataset",1:K,sep="") # naming studies
	} else {
    colnames(stat)<-colnames(P)<-names(data) 
	}
  if (!is.null(nperm)) {
    colnames(P.perm)<-paste("dataset",1:K,sep="") # naming studies
  }

  all.res<-list(stat=stat,p=P,bp=P.perm)
  attr(all.res$p,"nperstudy")<-N
  attr(all.res$p,"data.type")<- data.type
  attr(all.res$p,"individual.analysis")<-ind.method
  attr(all.res$p,"response.type")<- resp.type
  attr(all.res$p,"alter.hypothesis")<- tail
  return(all.res)
} #end of association with continuous phenotype
    
 if(resp.type == c("survival") )  {
  P<-stat<-P.perm<-NULL #prespecify the output matrix
  N<- NULL
  for (k in 1:K) {  
   y<- data[[k]] 
   r<- response.list[[k]]    
   ind.res<- get.p.logrank(y,r[,1], r[,2],asymptotic= asymptotic,
                         nperm=nperm,tail=tail)     
   P<- cbind(P,(ind.res$obs)[,2]) #p value
   stat<-cbind(stat,(ind.res$obs[,1])) # observed stats
   if (!is.null(nperm)) {
     P.perm<-cbind(P.perm,c(ind.res$perm)) #p value from permutation 
   }
   cat("dataset",k,"is done\n")
	 N<-c(N,ncol(y))	    
} # end of K study for loop

  if(is.null(names(data))) {
    colnames(stat)<-colnames(P)<-paste("dataset",1:K,sep="") #naming studies
	} else {
    colnames(stat)<-colnames(P)<-names(data)
	}
  if (!is.null(nperm)) {
    colnames(P.perm)<-paste("dataset",1:K,sep="") # naming studies
  }
  all.res<-list(stat=stat,p=P,bp=P.perm)
  attr(all.res$p,"nperstudy")<-N
  attr(all.res$p,"data.type")<- data.type
  attr(all.res$p,"individual.analysis")<-ind.method
  attr(all.res$p,"response.type")<- resp.type
  attr(all.res$p,"alter.hypothesis")<- tail
  return(all.res)    
} #end of association with survival 

} # end of Indi.DE.Analysis function 

## Output 
## twoclass: pvalue, log2FC and its SE (if available) 
## multiclass: pvalue 
## continuous, survival: stat, p, bp (perm.p)
       

get.sample.label.number<-function(lbl) {
  N<-length(lbl) #sample size per study 
  n<- table(lbl) #sample size per label per study
  return(list(N=N,n=n))
}

get.limma<-function(y,l,c,name,ANOVA,tail, paired){
## y: intensity matrix, l: group label, c: clinical data, name: group name
  if(!ANOVA) {
    if(is.null(c)) {
      if(paired) {
         subject <- as.factor(rep(1:(length(l)/2),2))
         design <-model.matrix(~ l + subject) 
      } else{
         design <-model.matrix(~ l)  # design matrix
      } 
    } else {
       if(paired) {
          subject <- as.factor(rep(1:(length(l)/2),2))
          lc <- data.frame(l=l,c=c,s=subject)
          design <-model.matrix(~ l + c + s,data=lc) 
       } else{
          lc <- data.frame(l=l,c=c)
          design <-model.matrix(~ l + c,data=lc) # design matrix
       }    
    }  
    #log2FC <- rowMeans(y[,l==name[1]])-rowMeans(y[,l==name[2]]) 
    fit <-lmFit(y, design)
    ebFit<-eBayes(fit)
    out.table <- topTable(ebFit,coef=2, number=Inf, sort.by='none')
    log2FC <- out.table$logFC
    lfcSE <- sqrt(ebFit$s2.post) * fit$stdev.unscaled[,2]
    #stat <- out.table$t
    p <- as.numeric(out.table$P.Value)
    dir <- sign(log2FC)
    if (tail=="high") { 
    	     pvalue <- mapply(function(x,y) if(x==1){
    	     	                             return(y/2) } else{
    	     	                             return(1-y/2)	
    	     	                         }, x=dir, y=p)
    	 }
    if (tail=="low") { 
    	 pvalue <- mapply(function(x,y) if(x==1){
    	     	                          return(1-y/2) } else{
    	     	                          return(y/2)	
    	     	                     }, x=dir, y=p)
     }    	 
     
    if (tail=="abs")  pvalue <- p  
    
    #pvalue <-out.table$P.Value
    limma.out <- list(log2FC,lfcSE,pvalue)
    names(limma.out) <- c('log2FC','lfcSE','pvalue')
  } else {
    design <-model.matrix(~ -1+l)  # design matrix
    colnames(design)[1:nlevels(l)] = make.names(levels(l))
    combination <- apply(combn(make.names(levels(l)),2),2,function(x) {
                                paste(x,collapse='-')
                              } )
    cont.matrix<-makeContrasts(contrasts=combination,levels = colnames(design))
    fit <-lmFit(y, design)
    fit2 <-contrasts.fit(fit, contrasts=cont.matrix)
    ebFit<-eBayes(fit2)
    pvalue <- ebFit$F.p.value
    limma.out <- list(pvalue)
    names(limma.out) <- c('pvalue')
  }
  return(limma.out)
}

get.sam<-function(y,l,name, seed, ANOVA, tail, paired) {
## y: intensity matrix, l: group label, name: group name 
  options(verbose=F)
  if(!ANOVA) {
    #log2FC <- rowMeans(y[,l==name[1]])-rowMeans(y[,l==name[2]]) 
    ## prepare the sam data
    sam.data<-list(x=y,y=l, geneid= 1:nrow(y), genenames=rownames(y),logged2=T)
    ## run the sam permutation test 
    samr.obj<-samr(sam.data,  resp.type=ifelse(paired,
                    "Two class paired","Two class unpaired"), nperms=100, 
                    random.seed=seed)
    #samr.obj<-samr(sam.data,  resp.type="Two class unpaired", nperms=100, 
    #               random.seed=seed)
    p<- as.numeric(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
    log2FC <- as.numeric(log2(samr.obj$foldchange))
    dir <- sign(log2FC)
    if (tail=="high") { 
    	     pvalue <- mapply(function(x,y) if(x==1){
    	     	                             return(y/2) } else{
    	     	                             return(1-y/2)	
    	     	                         }, x=dir, y=p)
    	 }
    if (tail=="low") { 
    	 pvalue <- mapply(function(x,y) if(x==1){
    	     	                          return(1-y/2) } else{
    	     	                          return(y/2)	
    	     	                     }, x=dir, y=p)
     }    	 
     
    if (tail=="abs")  pvalue <- p  
    
    sam.out <- list(log2FC,pvalue)
    names(sam.out) <- c('log2FC','pvalue')
  } else {
    sam.data<-list(x=y,y=l, geneid= 1:nrow(y), genenames=rownames(y),logged2=T)
    samr.obj<-samr(sam.data,resp.type="Multiclass", nperms=100,random.seed=seed)
    pvalue<- as.numeric(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
    sam.out <- list(pvalue)
    names(sam.out) <- c('pvalue')
  }
  options(verbose=T)
  return(sam.out)
}

get.limmaVoom <- function(y,l,c,name, ANOVA, tail, paired) {
## y: count matrix, l: group label, c: clinical data, name: group name 
  dge <- DGEList(counts=y)
  dge <- calcNormFactors(dge)  
  if(!ANOVA){
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject) 
      } else{
        design <-model.matrix(~ l)  # design matrix
      }  
    } else {
       if(paired) {
         subject <- as.factor(rep(1:(length(l)/2),2))
         lc <- data.frame(l=l,c=c,s=subject)
         design <-model.matrix(~ l + c + s,data=lc) 
       }else{
         lc <- data.frame(l=l,c=c)
         design <-model.matrix(~ l + c,data=lc) # design matrix
      }    
    }  
  v <- voom(dge,design,plot=FALSE,normalize="quantile") # voom normalization
  fit <-lmFit(v, design)
  ebFit<-eBayes(fit)
  out.table <- topTable(ebFit,coef=2, number=Inf, sort.by='none')
  log2FC <- out.table$logFC
  lfcSE <- sqrt(ebFit$var.prior)
  #stat <- out.table$t
    p <- as.numeric(out.table$P.Value)
    dir <- sign(log2FC)
    if (tail=="high") { 
    	  pvalue <- mapply(function(x,y) if(x==1){
    	     	                             return(y/2) } else{
    	     	                             return(1-y/2)	
    	     	                         }, x=dir, y=p)
    	 }
    if (tail=="low") { 
    	 pvalue <- mapply(function(x,y) if(x==1){
    	     	                          return(1-y/2) } else{
    	     	                          return(y/2)	
    	     	                     }, x=dir, y=p)
     }    	 
     
    if (tail=="abs")  pvalue <- p  
  
  limmaVoom.out <- list(log2FC,lfcSE,pvalue)
  names(limmaVoom.out) <- c('log2FC','lfcSE','pvalue')
  } else {
   design <-model.matrix(~ -1+l) 
   colnames(design)[1:nlevels(l)] = make.names(levels(l))
   combination <- apply(combn(make.names(levels(l)),2),2,function(x) {
                                paste(x,collapse='-')
                          } )
   cont.matrix<-makeContrasts(contrasts=combination,levels = colnames(design) )
   v <- voom(dge,design,plot=FALSE,normalize="quantile") # voom normalization
   fit <-lmFit(v, design)
   fit2 <-contrasts.fit(fit, contrasts=cont.matrix)
   ebFit<-eBayes(fit2)
   pvalue <- ebFit$F.p.value
   limmaVoom.out <- list(pvalue)
   names(limmaVoom.out) <- c('pvalue')
  }
  return(limmaVoom.out)
}

get.edgeR <- function(y,l,c,name, ANOVA, tail, paired) {
## y: count matrix, l: group label, c: clinical data, name: group name 
  dat <- DGEList(counts=y)
  dat <- calcNormFactors(dat)
  dat=estimateGLMCommonDisp(dat) #estimate common dispersion
  dat=estimateGLMTrendedDisp(dat) #estiamte trended dispersion
  dat=estimateGLMTagwiseDisp(dat) #estimate tagwise dispersion
  #dispersion <- dat$tagwise.dispersion
  
  if(!ANOVA){
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject) 
      } else{
        design <-model.matrix(~ l)  # design matrix
      }  
    } else {
       if(paired) {
         subject <- as.factor(rep(1:(length(l)/2),2))
         lc <- data.frame(l=l,c=c,s=subject)
         design <-model.matrix(~ l + c + s,data=lc) 
       }else{
         lc <- data.frame(l=l,c=c)
         design <-model.matrix(~ l + c,data=lc) # design matrix
       }    
    }      
    fit <- glmFit(dat, design)
    lrt <- glmLRT(fit, coef=2)
    out.table <-  lrt$table
    log2FC <- out.table$logFC
    
    p <- as.numeric(out.table$PValue)
    dir <- sign(log2FC)
    if (tail=="high") { 
    	 pvalue <- mapply(function(x,y) if(x==1){
    	     	                             return(y/2) } else{
    	     	                             return(1-y/2)	
    	     	                         }, x=dir, y=p)
    	 }
    if (tail=="low") { 
    	 pvalue <- mapply(function(x,y) if(x==1){
    	     	                          return(1-y/2) } else{
    	     	                          return(y/2)	
    	     	                     }, x=dir, y=p)
     }    	 
     
    if (tail=="abs")  pvalue <- p  
    
    edgeR.out <- list(log2FC,pvalue)
    names( edgeR.out) <- c('log2FC','pvalue')
   } else {
    if(is.null(c)) {
      design <-model.matrix(~ l)  # design matrix
    } else{
      lc <- data.frame(l=l,c=c)
      design <-model.matrix(~ l + c,data=lc) # design matrix
    }      
    fit <- glmFit(dat, design)
    lrt <- glmLRT(fit, coef=2:(nlevels(l)))
    out.table <-  lrt$table
    #log2FC <- out.table$logFC
    pvalue <- out.table$PValue
    edgeR.out <- list(pvalue)
    names(edgeR.out) <- c('pvalue') 
  }
  return(edgeR.out)
}

get.DESeq2 <- function(y,l,c,name, ANOVA, tail, paired) {
## y: count matrix, l: group label, c: clinical data, name: group name 
  #library(SummarizedExperiment)  
  if(!ANOVA){   
    if(is.null(c)) {
      if(paired) {
        subject <- as.factor(rep(1:(length(l)/2),2))
        design <-model.matrix(~ l + subject)  # design matrix
        colData <- data.frame(l=l, s=subject)
        colnames(colData) <-colnames(design)[-1]         
      }else{
        design <-model.matrix(~ l)  # design matrix
        colData <- data.frame(l=l)
        colnames(colData) <-colnames(design)[-1] 
      }  
    } else {
       if(paired) {
          subject <- as.factor(rep(1:(length(l)/2),2))
          lc <- data.frame(l=l,c=c,s=subject) 
          design <-model.matrix(~ l + c + s,data=lc)  # design matrix
          colData <- lc
          colnames(colData) <-colnames(design)[-1]           
       } else{
          lc <- data.frame(l=l,c=c) 
          design <-model.matrix(~ l + c,data=lc)  # design matrix
          colData <- lc
          colnames(colData) <-colnames(design)[-1] 
       } 
    }
      
  ddsMat <- DESeqDataSetFromMatrix(countData = y,
                                    colData = colData,
                                    design = as.formula(
                        paste(" ~ ",paste(colnames(colData), collapse=" + ") ) 
                        )  
                      )
  ddsMat <- DESeq(ddsMat)
  res <- results(ddsMat,contrast=c(colnames(colData)[1],levels(l)[2],
                                   levels(l)[1]) )
    log2FC <- as.numeric(res$log2FoldChange)
    lfcSE <- as.numeric(res$lfcSE)
    
    p <- as.numeric(res$pvalue)
    dir <- sign(log2FC)
    if (tail=="high") { 
    	  pvalue <- mapply(function(x,y) if(x==1){
    	     	                             return(y/2) } else{
    	     	                             return(1-y/2)	
    	     	                         }, x=dir, y=p)
    	 }
    if (tail=="low") { 
    	 pvalue <- mapply(function(x,y) if(x==1){
    	     	                          return(1-y/2) } else{
    	     	                          return(y/2)	
    	     	                     }, x=dir, y=p)
     }    	 
     
    if (tail=="abs")  pvalue <- p  
    
    DESeq2.out <- list(log2FC,lfcSE,pvalue)
    names(DESeq2.out) <- c('log2FC','lfcSE','pvalue')  
  } else {     
    if(is.null(c)) {
       design <-model.matrix(~ l) 
       colData <- data.frame(group=l)
    } else{
       lc <- data.frame(l=l,c=c)
       design <-model.matrix(~ l + c,data=lc)
       colData <- lc
       colnames(colData) <-colnames(design)[-1] 
    }
    
  ddsMat <- DESeqDataSetFromMatrix(countData = y,
                                     colData = colData,
                                     design = as.formula(
                        paste(" ~ ",paste(colnames(colData), collapse=" + ") )
                      )  
                    )    
  ddsMat <- DESeq(ddsMat,test='LRT',reduced= as.formula(
     paste(" ~ ",paste(colnames(colData)[-c(1:(nlevels(l)-1))], collapse=" + "))  
    )  
  )
    res <- results(ddsMat)
    pvalue <-  as.numeric(res$pvalue)
    DESeq2.out <- list(pvalue)
    names(DESeq2.out) <- c('pvalue')      
  }  
  return(DESeq2.out)
}

#--------calculate r statistic (pearson's correlation)for all genes-----#
#---note: l is continuous------------#
cal.pearsonr<-function(y,l) {
  stopifnot(length(l)==ncol(y), is.matrix(y))
  n<-length(l)
  num<-n*y%*%l-rowSums(y)*sum(l)
  den<-sqrt(n*rowSums(y^2)-rowSums(y)^2)*sqrt(n*sum(l^2)-sum(l)^2)
  r<-num/den
  r
}


#-------get p value for pearson correlation r using permutation-------#

get.p.pearsonr<-function(y,l,asymptotic, nperm,tail) {
## y: expression matrix , l: continuous outcome, tail: specify the Ha
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-q<-rep(NA,nrow(y))
    names(stat.obs)<-names(p)<-rownames(y)
    stat.obs[rnum]<-cal.pearsonr(y[rnum,],l) #observed stat
   if (asymptotic==F) {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
    stat.perm[rnum,]<-replicate(nperm,cal.pearsonr(y[rnum,],l)) #stat from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
   } else {
     n<-length(l)
     t<-stat.obs*sqrt((n-2)/(1-stat.obs^2)) #fisher's transformation
     if (tail=="low") p[rnum]<-pt(t,n-2)
     if (tail=="high") p[rnum]<-1-pt(t,n-2)
     if (tail=="abs") p[rnum]<-2*(pmin(pt(t,n-2),1-pt(t,n-2)))
     res<-list(obs=cbind(stat.obs,p))
   }
  return(res)
}

#--------calculate r statistic (spearman's correlation)for all genes------#
#---note: l is continuous------------#
cal.spearmanr<-function(y,l) {
    stopifnot(length(l)==ncol(y), is.matrix(y))
    n<-length(l)
    y<-t(apply(y,1,rank))
    l<-rank(l)
    num<-n*y%*%l-rowSums(y)*sum(l)
    den<-sqrt(n*rowSums(y^2)-rowSums(y)^2)*sqrt(n*sum(l^2)-sum(l)^2)
    r<-num/den
    r
}

#-------get p value for spearman correlation r using permutation----------#
get.p.spearmanr<-function (y,l,asymptotic, nperm,tail) {
## y: expression matrix , l: continuous outcome, tail: specify the Ha
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-q<-rep(NA,nrow(y))
    names(stat.obs)<-names(p)<-rownames(y)

    stat.obs[rnum]<-cal.spearmanr(y[rnum,],l)#observed stat

   if (asymptotic==F)
   {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)

    stat.perm[rnum,]<-replicate(nperm,cal.spearmanr(y[rnum,],l))# stat from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
   } else {
     n<-length(l)
     t<-stat.obs*sqrt((n-2)/(1-stat.obs^2)) #fisher's transformation
     if (tail=="low") p[rnum]<-pt(t[rnum],n-2)
     if (tail=="high") p[rnum]<-1-pt(t[rnum],n-2)
     if (tail=="abs") p[rnum]<-2*(pmin(pt(t[rnum],n-2),1-pt(t[rnum],n-2)))
     res<-list(obs=cbind(stat.obs,p))
   }
  return(res)
}

#-------get logrank z statistic-----------------#
cal.logrank<-function(y,time,event)
{
  get.stat<-function(x,time,event) {  
    stat<-summary(coxph(Surv(time,event)~x),method="breslow")$sctest[1]
    stat
  }
  z<-apply(y,1,get.stat,time=time,event=event)
  z
}
cal.p.logrank<-function(y,time,event)
{
  get.p<-function(x,time,event) {  
    p<-summary(coxph(Surv(time,event)~x,method="breslow"))$sctest[3]
    p
  }
  p<-apply(y,1,get.p,time=time,event=event)
  p
}

#-------get p value for logrank z using permutation-----------------#
get.p.logrank<-function(y,time,event,asymptotic,nperm,tail)
{
  rnum<-which(apply(y,1,function(x) !any(is.na(x))))
  stat.obs<-p<-pp<-rep(NA,nrow(y))
  
  stat.obs[rnum]<-cal.logrank(y[rnum,],time=time,event=event)
  
  names(stat.obs)<-names(p)<-rownames(y)
  if(asymptotic== F)
  {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
    stat.perm[rnum,]<-replicate(nperm,cal.logrank(y[rnum,],time=time,
                                                  event=event))
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail)
    #q[rnum]<-p.adjust(p[rnum],method="BH")    
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
  } else {
    p[rnum]<-cal.p.logrank(y[rnum,],time=time,event=event)
    if (tail=="low")  pp[rnum]<-p
    if (tail=="high") pp[rnum]<-1-p
    if (tail=="abs")  pp[rnum]<-2*(pmin(p,1-p))
    res<-list(obs=cbind(stat.obs,p=pp))
  }
  res
}



