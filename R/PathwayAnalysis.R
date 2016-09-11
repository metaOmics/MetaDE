##' Main Function for pathway analysis
##' The \code{PathAnalysis} is a function to perform pathway analysis 
##' (a.k.a. gene set enrichment test) for functional annotation of a candidate
##' gene list or an ordered gene result from Meta DE analysis output.
##' @title Main Function for pathway analysis.  
##' @param meta.p is a vector of meta-analysis p-value.
##' @param pathway is a vector of pathway databases used for functional 
##' analysis, see data(pathways) for more details.
##' @param enrichment is the method used for pathway analysis, must be one of
##' "KS" and "Fisher's exact".
##' @param p.cut is the p-value cutoff to select the DE genes, option for 
##' Fisher's exact method only. 
##' @param DEgene.number is the top number of DE genes, option for 
##' Fisher's exact method only. 
##' @param size.min is the minimum pathway size to be included in the 
##' functional analysis. 
##' @param size.max is the maximum pathway size to be included in the 
##' functional analysis. 

##' @return a data frame with columns: \cr
##' \itemize{
##' \item{pvalue}{ the p-value from pathway analysis for each pathway}
##' \item{qvalue}{ the q-value from pathway analysis for each pathway}
##' \item{OddsRatio}{ optional, the odds ratio from Fisher's exact test method}
##' \item{logOR}{ optional, the log odds ratio from Fisher's exact test method}
##' \item{DEgenes}{ optional, the set of DE genes in each pathway}
##' }

##' @export
##' @examples
##' \dontrun{
##' meta.p  <- meta.res$meta.analysis$pval
##' ks.result <- PathAnalysis(meta.p = meta.p, enrichment = "KS")
##' fisher.result <- PathAnalysis(meta.p = meta.p, enrichment = "Fisher's exact")
##' }

PathAnalysis <- function (meta.p = NULL,pathway = c(Biocarta.genesets,
                           GOBP.genesets,GOCC.genesets,GOMF.genesets,
                           KEGG.genesets,Reactome.genesets),
                           enrichment = c("KS","Fisher's exact"),
                           p.cut=NULL,DEgene.number = 200,
                           size.min = 15, size.max = 500)
{ 
  if(is.vector(meta.p)) {
  	all.genes <- names(meta.p)
  	meta.p <- matrix(meta.p)
  	rownames(meta.p) <- all.genes
  	colnames(meta.p)[1] <- "pvalue"
  }	
  		
  #library(MetaPath)
  data(pathways)	
  enrichment = match.arg(enrichment)
  pathway = pathway[which(sapply(pathway,length) >= size.min & 
                             sapply(pathway,length) <= size.max)]
  gene.in.DB = unique(unlist(pathway))
  set.name = names(pathway)
  gene.in.array = rownames(meta.p)
  gene.common = intersect(gene.in.array, gene.in.DB)
  DB.matrix = matrix(0, length(set.name), length(gene.common))
    rownames(DB.matrix) = set.name
    colnames(DB.matrix) = gene.common
    colnames(DB.matrix) = toupper(colnames(DB.matrix))
  for (t1 in 1:length(set.name)) {
    gene = toupper(intersect(pathway[[t1]], gene.common))
    DB.matrix[set.name[t1], gene] = 1
  }
  rm(pathway)
  if (enrichment == "KS"){
    gene.name.sort = names(sort(meta.p[,1],decreasing = F))
    gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
    gene.name.sort = toupper(gene.name.sort)
    set2allgenes.mtx = DB.matrix
    order.mtx.1 = (set2allgenes.mtx[, gene.name.sort])
    order.mtx.0 = (1 - order.mtx.1)
    n_hit = rowSums(order.mtx.1)
    n_miss = rowSums(order.mtx.0)
    n_genes = ncol(order.mtx.1)
    nn = (n_hit*n_miss)/n_genes
    order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
    order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
    order.mtx = order.mtx.0 + order.mtx.1
    ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                           max))
    pvalue.0 = matrix(exp(-2*nn*(ES.0[,1]^2)), ncol = 1)
    rownames(pvalue.0) = rownames(ES.0)
    pvalue.set.0 = pvalue.0 
    qvalue.set.0 = as.matrix(p.adjust(pvalue.0,"BH"),ncol = 1)
    rownames(qvalue.set.0) = rownames(pvalue.0)
    out <- data.frame(pvalue=pvalue.set.0,qvalue=qvalue.set.0)
    return(out)
    #return(list(pvalue.meta = pvalue.set.0, qvalue.meta = qvalue.set.0))
   }
  else if (enrichment == "Fisher's exact")
   {
   	if (is.null(p.cut)) { 
      gene.name.sort = names(sort(meta.p[,1],decreasing = F))
      genes.in.study = gene.name.sort
      gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
      gene.name.sort = toupper(gene.name.sort)
      DEgene = gene.name.sort[1:DEgene.number]
     } else {
      gene.name.sort = names(sort(meta.p[,1],decreasing = F))
      genes.in.study = gene.name.sort
      gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
      gene.name.sort = toupper(gene.name.sort)     	
      gene.select = names(which(meta.p[,1] < p.cut))
      DEgene = toupper(gene.select[gene.select%in%gene.common])
     }      
     
    DEset = logOR = OR = pvalue.0 = matrix(NA,nrow = nrow(DB.matrix),ncol = 1)
    rownames(pvalue.0) = rownames(DB.matrix)
    for (i in 1:nrow(DB.matrix)){
      count_table<-matrix(0,2,2)
      p_value <-NA
      ####in the gene list and in the pathway
      count_table[1,1]<-sum(DEgene %in% colnames(DB.matrix[,DB.matrix[i,]==1]))
      ####in the gene list but not in the pathway
      count_table[1,2]<-length(DEgene)-count_table[1,1]
      ####not in the gene list but in the pathway
      count_table[2,1]<-sum(genes.in.study%in% 
                               colnames(DB.matrix[,DB.matrix[i,]==1]))
      ####not in the gene list and not in the pathway
      count_table[2,2]<-length(genes.in.study)-count_table[2,1]       
     if(length(count_table)==4){
       pvalue.0[i,1] <- fisher.test(count_table, alternative="greater")$p}
       OR[i,1] <- (count_table[1,1] * count_table[2,2])/
                       (count_table[1,2] * count_table[2,1])
       logOR[i,1] <- ifelse(OR[i,1]==0,0,log(OR[i,1]))
       DEset[i,1] <- paste(DEgene[DEgene %in% 
                       colnames(DB.matrix[,DB.matrix[i,]==1])], collapse= ",")
    }
    OR.0 = OR[pvalue.0[,1]!=1,,drop=FALSE]
    logOR.0 = logOR[pvalue.0[,1]!=1,,drop=FALSE]
    DEset.0 = DEset[pvalue.0[,1]!=1,,drop=FALSE]
    pvalue.0 = pvalue.0[pvalue.0[,1]!=1,,drop=FALSE]
    qvalue.0 = matrix(p.adjust(pvalue.0[,1], "BH"),ncol = 1)
    rownames(qvalue.0) = rownames(pvalue.0)
    pvalue.set.0 = pvalue.0 
    qvalue.set.0 = as.matrix(p.adjust(pvalue.0,"BH"),ncol = 1)
    rownames(qvalue.set.0) = rownames(pvalue.0)
        
    out <- data.frame(pvalue=pvalue.set.0,qvalue=qvalue.set.0, 
                      OddsRatio=OR.0, logOR=logOR.0, DEgenes=DEset.0)
    return(out)
    #return(list(pvalue.meta = pvalue.set.0, qvalue.meta = qvalue.set.0))
   }
}

## example  
#xx = MAPE_G(meta.p = meta.res.p$meta.analysis$pval,enrichment = "KS")
#yy = MAPE_G(meta.p = meta.res.p$meta.analysis$pval,enrichment = "Fisher's exact")


