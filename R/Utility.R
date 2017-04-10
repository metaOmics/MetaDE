## Utility functions for plotting and other non stats purposes (Internal)
## Author: Tianzhou Ma
## Institution: University of pittsburgh
## Date: 08/29/2016

#----------------------------------------------# 
# Matrix manipulation methods 
#----------------------------------------------# 
# Flip matrix (upside-down) 
flip.matrix <- function(x) { 
    mirror.matrix(rotate180.matrix(x)) 
}   

# Mirror matrix (left-right) 
mirror.matrix <- function(x) { 
    xx <- as.data.frame(x); 
    xx <- rev(xx); 
    xx <- as.matrix(xx); 
    xx; 
} 

# Rotate matrix 180 clockworks 
rotate180.matrix <- function(x) { 
    xx <- rev(x); 
    dim(xx) <- dim(x); 
    xx; 
} 

#-----------------------------------------------------#
# generate nPr    with repetition                     #
# for generating all possible weights                 #
#-----------------------------------------------------#  
permut<-function (n, r) { 
    v<-1:n
    sub <- function(n, r, v) {
        if (r == 1) matrix(v, n, 1)
        else if (n == 1) matrix(v, 1, r)
        else {
            inner <- Recall(n, r - 1, v)
            cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
            ncol = ncol(inner), nrow = nrow(inner) * n, byrow = TRUE))
       }
    }
   sub(n, r, v[1:n])
}
#-----------------------------------------------------#
#        wt=possible weights for each datasets
#        n=# of datasets                              #
#-----------------------------------------------------#
gen.weights<-function(wt,n) {
    comb<-permut(n=length(wt),r=n)
    weight<-matrix(wt[comb],ncol=n)
    return(weight[-1,])
}

#--------------------------------------------------------------#
# plot a matrix of gene expression data with rows are genes
# columns are samples
# n: # of samples in each study
# ni: # of samples in each class of each study
#--------------------------------------------------------------#
plot.matrix<-function(mat,color="GR") {
    n<-attr(mat,"n")
    ni<-attr(mat,"ni")
    label<-attr(mat,"label")
    dataset.name <- attr(mat,"dataset.name")
    group.name <- attr(mat,"group.name")
    ref.group <- attr(mat,"ref.group")
    nc<-ncol(mat)
    nr<-nrow(mat)
    cexCol1<-1/log(nc)
    cexCol2<-1/log10(nc)
    cexRow<-1/log(nr,4)
    K<-length(n)
    #mycol<- c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", 
    #"#238443", "#006837", "#004529")
    colfun<-colorRampPalette(c("green","black","red"))
    if(color=="BY") colfun<-colorRampPalette(c("blue","black","yellow"))
    mycol<-colfun(16)
    mval<-min(max(abs(mat),na.rm=T),3)
    xcut<-seq(-mval,mval,len=length(mycol)-1)
    xcut<-c(-100,xcut,100)
    m <- matrix(1:2, 2, 1)
    nf<-layout(m, heights=c(6, 1))
    par(mar=c(2,3,2,5))
    image(x=1:nc,y=1:nr,z=t(flip.matrix(mat)),col=mycol,axes=FALSE,xlab = "", 
          ylab = "", breaks=xcut)
    #axis(3,1:nc,labels=label,las=2,line=-0.5,tick=0,cex.axis=cexCol1)
    axis(3,cumsum(ni)-ni/2+0.5,labels=rep(group.name,K),las=1,line = -0.5,
         tick=0,cex.axis=cexCol2)
    axis(4,nr:1,labels=(row.names(mat)), las = 2, line = -0.5, tick = 0, 
         cex.axis = cexRow)
    axis(1,cumsum(n)-n/2+0.5,labels=dataset.name,las=1,line = -1,
         tick=0,cex.axis=cexCol2)
    #---distinguish studies----#
    abline(v=cumsum(n)+0.5,lwd=2,col="white")
    #---distinguish classes----#
    if (is.null(ncol(label))) abline(v=cumsum(ni)+0.5,lwd=2,col="white",lty=2)
    #---------if AW method we add category information on the plot---------#
    if (!is.null(attr(mat,'category'))) {
        cat<-attr(mat,'category')
        at.line <-cumsum(table(cat))+0.5
        sum.pos <- cumsum(table(cat))
        at <- rep(NA, length(sum.pos))
        for (i in 2:length(sum.pos)){
        	at[i] <- sum.pos[i-1] + (sum.pos[i] - sum.pos[i-1])/2 + 0.5
          }
        at[1] <- sum.pos[1] + (sum.pos[1] - 0)/2 + 0.5
        axis(2,at-0.5,labels=rev(unique(cat)),tick = 0, las=1,
             cex.axis = cexRow+0.2)
        abline(h=at.line,lwd=2,col="white")
    }
    #----add legend---------------#
    l<-length(xcut)
    image(1:(l-1),0.5,matrix(xcut[-1],nrow=l-1,ncol=1),col=mycol,breaks=xcut,
          axes=F,xlab="",ylab="")
    marcas<-(0:(l-1))+0.5
    axis(1,marcas,round(xcut,1),tick=0.5,cex.axis=cexCol2,line=-0.5)
}

#-----------------------------------------------------------#
# order genes for plotting obtained from all methods except AW  #
#-----------------------------------------------------------#
order.genes.simple<-function(dat) {
    r<-hclust(dist(dat))$order
    plot.genes<-dat[r,]
    return(plot.genes)
}

#---------------------------------------------------------------------------#
# output genes according to the order of categories defined from AW.weight  #
#---------------------------------------------------------------------------#
order.genes.AW<-function(dat,AW.weight) {
        K<-ncol(AW.weight)
        #wt.group<-gen.weights(c(0,1),K) # generate all possible weights
		wt.group<-do.call(expand.grid, rep(list(c(0, 1)), K))[-1,]
        wt.group<-wt.group[order(apply(wt.group,1,sum)),]
        row.names(wt.group)<-nrow(wt.group):1 

        ng<-nrow(dat) # number of significant genes
        group<-rep(NA,ng)

        #---------order the OW categories nicely ------------------#
        for (i in 1:ng) {
            for (j in 1:nrow(wt.group)) {
                if (sum(AW.weight[i,]==wt.group[j,])==K) {
                    group[i]<-row.names(wt.group)[j]
                    next
                }
            }
        }
        #perform hierarchical clustering in each weight group
        plot.genes.ordered<-NULL
        for (i in sort(as.numeric(names(table(group))))) {
            x<-subset(dat,group==i)
            if (nrow(x)>2) {
                newx<-x[hclust(dist(x))$order,]
            } else newx<-x
            plot.genes.ordered<-rbind(plot.genes.ordered,newx)
        }
    attr(plot.genes.ordered,"category")<-
      apply(wt.group,1,paste,collapse=',')[sort(group)]
    return(plot.genes.ordered)
}

#-----------------------------------------------------------------------------#
# check gene names
#-----------------------------------------------------------------------------#
check.exp<-function(x)
{
  if (is.null(row.names(x[[1]]))) {
   K<-length(x)
   ng<-nrow(x[[1]])
   for (k in 1:K) {
     row.names(x[[k]])<-paste("gene",1:ng)
   }
 }
 return(x)
}

#------------------------------------------------------------------------------#
# check dimensions and size of argument
#------------------------------------------------------------------------------#
check.dim<-function(x,y,ind.method,meta.method,paired){
	K<-length(x)
	nperstudy<-sapply(x,function(y) ncol(y))
	nlabels<-sapply(y,function(z) length(z))						
	if(sum(nperstudy==nlabels)!=K) {
    stop(cat("The number of samples does not match with the dimension of 
             labels in study(s)",paste((1:K)[nperstudy!=nlabels],"",
                                       collapse=","),"!"))
	}
	if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
		if(length(ind.method)!=K)stop(paste('Argument "ind.method" should be 
                                        a character vector of size',K))
	}
	if(("REM"%in%meta.method|"FEM"%in%meta.method)&(length(paired)!=K)) {
    stop(paste('Argument "paired" should be a logical vecter of size',K))
	}
}

#-----------------------------------------------------------------------------#
#   check if parametric is ok                                                 #
#-----------------------------------------------------------------------------#
check.parametric<-function(meta.method,parametric)
{
  if (parametric==TRUE&sum(meta.method%in%c("SR","PR","rankProd",
                                        "Fisher.OC","minMCC"))>0) {
    stop(paste("There is no parametric result for",meta.method))
  }
}

#-----------------------------------------------------------------------------#
#   check tail                                         #
#-----------------------------------------------------------------------------#
check.tail<-function(meta.method,tail)
{
  if (tail=='abs'&length(grep('.OC',meta.method))>0) {
    stop(paste("If you chose",meta.method,",then you should specify the 'tail' 
               to be either 'high' or 'low'"))
  }
}

#-----------------------------------------------------------------------------#
#   check individual methods and response                                       
#-----------------------------------------------------------------------------#
check.indmethod<-function(y, resp.type, data.type,ind.method, 
                             select.group,tail,paired) {  
  K<-length(y)
for(k in 1:K) {
  if(data.type == "continuous") {
   if(!(ind.method[k] %in%c('limma','sam',"pearsonr","spearmanr",'logrank'))) {
    stop ("Incorrect method for microarray or RNAseq FPKM")
   }
  } else if (data.type =='discrete') {
     if(!(ind.method[k] %in%c('edgeR','DESeq2',"spearmanr",'limmaVoom') ) ) {
       stop ("Incorrect method for RNAseq count")
     } 
    }   
  if(resp.type %in% c("twoclass", "multiclass")) {
    if(!(ind.method[k] %in%c("limma","sam","limmaVoom","edgeR","DESeq2") ) ) {
      stop ("Incorrect method for the response")
    }  
  }
  if(resp.type =="twoclass") {
    if(nlevels(as.factor(y[[k]]))!=2 && length(select.group)!=2 ) {
      stop (cat(resp.type," requires two levels ") )
    }
    if (is.null(paired)) {
      stop(paste("you need to specify the logical value for 'paired' 
               for study", k))
    }
  }     
  if(resp.type=="multiclass") {
    if(nlevels(as.factor(y[[k]])) <=2) {
      stop (cat(resp.type," requires at least three levels") )
    }
    if(tail!="abs") {
      stop (cat(resp.type," cannot perform one-sided test") )
    }
  }  
  if(resp.type %in% c("continuous") ) {
    if(!(ind.method[k] %in%c("pearsonr","spearmanr") ) ) { 
      stop ("Incorrect method for the response") 
    }     
    if(!is.numeric(y[[k]])) {
      stop (cat(resp.type," response has to be numeric") )
    }
  }
  else if (resp.type %in% c("survival") ) {
    if(!(ind.method[k] %in%c("logrank") ) ) {
      stop ("Incorrect method for the response")
    }
    if(is.null(y[[k]][,1])) {
      stop( cat("dataset",k, "is missing the event time")) 
    }
    if(is.null(y[[k]][,2])) {
      stop( cat("dataset",k, "is missing the censoring status"))
    } 
    
    if(!is.numeric(y[[k]][,1])) {
      stop( cat("dataset", k, ": Survival time has to be numeric")) 
    }
    
    if(!(y[[k]][,2] %in% c(0,1)| y[[k]][,2] %in% c('0','1') )) {
      stop( cat("dataset", k, ": Censoring status has to be either 0 or 1")) 
    }
      
  }
 }  
}

#----------------------------------------------------------------------------#
# check meta methods
#----------------------------------------------------------------------------#

check.metamethod<- function(x,y, resp.type,ind.method,meta.method,rth=NULL,
                            REM.type=NULL,paired=NULL) {
  ## x is the data, y is the response
  K<-length(x)
  cat("Please make sure the following is correct:\n")
  cat("*You input",K,"studies\n")
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
    cat("*You selected",ind.method,"for your",K,"studies respectively\n")
  }
  if (sum(paired)>0) cat("*Some of the studies are paired design\n")
  if (is.null(paired)) cat("*They are not paired design\n")
  cat("*",meta.method, "was chosen to combine the",K,"studies,respectively\n")
  if (length(meta.method)>1 & 
        sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))>0) {
    stop("Sorry, we currently do not allow multiple choices of meta.method 
         for 'FEM','REM','rankProd','minMCC'")
  }
  
  if (resp.type %in% c("continuous", "survival")) {
  	   if (meta.method %in% c("FEM","REM","minMCC","rankProd") ) {
         stop("The meta.methods of 'FEM','REM','rankProd','minMCC' are not 
              applicable to 'continuos','survival' response") 
  	   }
  } 
  
#  if(sum(meta.method%in%c("maxP.OC","minP.OC","Fisher.OC","roP.OC",
#                          "Stouffer.OC"))>0) {
#  	  if (tail == c("abs") ) {
#         stop("One-side corrected method ") 
#  	   }
#  }	
  	  
  if(sum(meta.method%in%c("FEM","REM","minMCC","rankProd"))<1) {
    for(i in 1:K){
      check.label.for.meta(y[[i]],method=ind.method[i], 
                           meta.method=meta.method,k=i)
    }
  }
  
  if (length(grep('roP',meta.method))!=0&is.null(rth)) {
    stop("You should specify rth=XXX, when you choose roP method")
  }
  if (length(grep('REM',meta.method))!=0&is.null(REM.type)) {
    stop("You should specify REM.type=XXX, when you choose REM method")
  }
  if (!is.null(rth)){
    if (length(grep('roP',meta.method))!=0&&length(x)<rth) {
      stop("rth shouldn't be larger than the number of datasets")
    }
  } 
  
  if (!is.null(paired)&length(paired)<K) {
    stop(paste("you need to specify a vector of logical value for 'paired' 
               for all",K,"studies"))
  }
  
  ANOVA <- (resp.type == "multiclass")
  if (length(grep('maxP',meta.method))==0&&ANOVA) {
    warning("maxP is suggested for ANOVA model")
  }
}

check.label.for.meta<-function(L,method,meta.method,k)
{
  if(!is.null(dim(L)))stop(cat("Please check whether the dimension in study",k,
                               " is matched to individual method",method,"?"))
 nL<-nlevels(as.factor(L))
 if ( nL<2) stop("All individual methods requires at least two levels")
 if (length(grep('minMCC',meta.method))!=0)
  {
   if (nL==2) warning("minMCC could test for two classes. But for better 
                      performance, try limma+maxP or sam+maxP")
   if (nL<2) stop(paste(meta.method,"minMCC method requires at least two 
                        levels"))
  }
 if (sum(table(L)<=1)==1) stop("<= one sample in the group, can not do the test 
                               or check labels")
 if (!is.null(method)&length(grep('minMCC',meta.method))!=0) 
    warning(paste("minMCC is a method that can not be combined with",method,". 
                  We'll perform minMCC only."))
 if (!is.null(method)&length(grep('rankProd',meta.method))!=0) 
    warning(paste("rankProd is a method that can not be combined with",method,".
                  We'll perform rankProd only."))

}


#------------------------------------------------------------------------------#
#perm.lab: function to permute the labels of disease status
# x: labels
# paired: a logical to specify whether the data is pair-designed or not
#------------------------------------------------------------------------------#
perm.lab<-function(x,paired=FALSE){
	New.x<-x	
	if(paired){
		templab<-rbinom(length(x)/2,1,0.5)
		index.d<-which(x==levels(x)[2])
		index.c<-which(x==levels(x)[1])
		dx<-x[index.d]
		cx<-x[index.c]
		New.dx<-dx
		New.cx<-cx
		New.dx[which(templab==1)]<-cx[which(templab==1)]
		New.cx[which(templab==1)]<-dx[which(templab==1)]
		New.x[index.d]<-New.dx
		New.x[index.c]<-New.cx
	}else{
		New.x<-sample(x)
	}
	return(New.x)
}

#=============================================================================#
# summary the DE number in a table
# pm: the p-value matrix
# p.cut: a numeric vector of p-values at which the DE numbers are counted 
# q.cut: a numeric vector of q-values at which the DE numbers are counted
# method: a vector of character string specifying the method
#-----------------------------------------------------------------------------#
count.DEnumber<-function(result,p.cut,q.cut){
	  if(class(result)=="MetaDE.pvalue"){
			pm<-cbind(result$ind.p,result$meta.analysis$pval) 
		}else if(class(result)=="MetaDE.ES"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else if(class(result)=="MetaDE.minMCC"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else{
			pm<-result
		}
        qm<-cbind(apply(pm,2,function(x)p.adjust(x,method="BH")))
        table.p<-matrix(NA,length(p.cut),ncol(pm))
        for(i in 1:length(p.cut)){
                table.p[i,]<-apply(pm,2,function(x)sum(x<=p.cut[i],na.rm=T))
        }
        table.q<-matrix(NA,length(q.cut),ncol(pm))
        for(i in 1:length(q.cut)){
                table.q[i,]<-apply(qm,2,function(x)sum(x<=q.cut[i],na.rm=T))
        }       
        rownames(table.p)<-paste("p=",p.cut,sep="")
        rownames(table.q)<-paste("FDR=",q.cut,sep="")
        colnames(table.p)<-colnames(table.q)<-colnames(pm)
        return(list(pval.table=table.p,FDR.table=table.q))
}

draw.DEnumber<-function(result,maxcut,mlty=NULL,mcol=NULL,mlwd=NULL,mpch=NULL,
                        FDR=TRUE){
		if(class(result)=="MetaDE.pvalue"){
			pm<-cbind(result$ind.p,result$meta.analysis$pval) 
		}else if(class(result)=="MetaDE.ES"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else if(class(result)=="MetaDE.minMCC"){
			pm<-cbind(result$meta.analysis$pval)
			colnames(pm)<-attr(result$meta.analysis,"meta.method")
		}else{
			pm<-result
		}
         method<-colnames(pm)
        if(FDR) pm<-cbind(apply(pm,2,function(x)p.adjust(x,method="BH")))
         maxp<-max(pm,na.rm=T)
         if (maxcut>maxp)
        {
          cat("Your maximum cut point exceeds the maximum of observed 
              p-value/FDR\n",
           "we will use",maxp,"as the maximum cut point\n")
           maxcut<-maxp
         }
	  ns<-ncol(pm)
        ymax<-max(apply(pm,2,function(x)sum(x<=maxcut,na.rm=T)))
        if(is.null(mlty))mlty=1:ns
        if(is.null(mcol))mcol=1:ns
        if(is.null(mlwd))mlwd=rep(2,ns)
        if(is.null(mpch))mpch=1:ns
        xlab0<-ifelse(FDR,"FDR cut-off","p-value cut-off")
       #----------get an optimal place to draw the symbols--------#
       get.c<-function(cut,pm)
       {
        s<-apply(pm,2,function(x,y) sum(x<=y,na.rm=T),y=cut)        
        return(sum(dist(cbind(cut,s))))
       }
       mycut<-as.matrix(seq(0,maxcut,length=20))
       dis<-apply(mycut,1,get.c,pm=pm)
       minx.pos<-mycut[which.max(dis)]
   
        plot(c(0,maxcut),c(1,ymax),type='n',xlab=xlab0,ylab="Significant tests")
        for(i in 1:ns){		
			y.pos<-sum(pm[,i]<=minx.pos,na.rm=T)
			if(y.pos==0){x.pos<-minx.pos
			}else{
            x.pos<-sort(pm[,i])[y.pos]}
			points(x.pos,y.pos,pch=mpch[i],col=mcol[i],lwd=3)
			lines(sort(pm[,i]),rank(sort(pm[,i]),ties.method="max"),lty=mlty[i],
            col=mcol[i],lwd=mlwd[i])
        
}
        legend("topleft",method,lty=mlty,lwd=mlwd,col=mcol,bty='n',pch=mpch)
}
