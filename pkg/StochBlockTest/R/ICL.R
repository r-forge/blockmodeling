# not exported
findActiveParam<-function(M, n, k, na.rm=TRUE){
  parByHB<-outer(k,k)
  parByHB<-array(parByHB,dim=c(dim(parByHB),dim(M)[3]))
  if(sum(n)!=dim(M)[1]) stop("Dimensions do not match!")
  kVec<-1:length(n)
  clu<-rep(kVec, n)
  for(r in 1:dim(M)[3]){
    for(i in kVec)for(j in kVec){
      parByHB[i, j,r]<-parByHB[i, j,r]*(var(as.vector(M[clu==i, clu==j,r]), na.rm = na.rm)>0)
    }
  }
  return(parByHB)
}



#' Function that computes integrated classification likelihood based on stochastic one-mode and linked block modeling. If \code{clu} is a list, the method for linked/multilevel networks is applied. The support for multirelational networks is not tested.
#'
#' @param M A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation.
#' @param clu A partition. Each unique value represents one cluster. If the nework is one-mode, than this should be a vector, else a list of vectors, one for each mode. Similarly, if units are comprised of several sets, clu should be the list containing one vector for each set.
#' @param weights The weights for each cell in the matrix/array. A matrix or an array with the same dimmensions as \code{M}. 
#' @param uWeights The weights for each unin. A vector with the length equal to the number of units (in all sets).
#' @param diagonal How should the diagonal values be treated. Possible values are:
#' \itemize{
#'   \item ignore - diagonal values are ignored 
#'   \item seperate - diagonal values are treated seperately
#'   \item same - diagonal values are treated the same as all other values
#' }
#' @param limits If \code{diagonal} is \code{"ignore"} or \code{"same"}, an array with dimmensions equal to:
#' \itemize{
#'   \item number of clusters (of all types)
#'   \item number of clusters (of all types)
#'   \item number of relations
#'   \item 2 - the first is lower limit and the second is upper limit
#' }
#' If \code{diagonal} is \code{"seperate"}, a list of two array. The first should be as described above, representing limits for off diagonal values. The second should be similar with only 3 dimensions, as one of the first two must be omitted.
#' 
#' @return The value of ICL
ICLStochBlock<-function(M, 
                       clu, 
                       weights=NULL,
                       uWeights=NULL,
                       diagonal = c("ignore","seperate","same"),
                       limitType=c("none","inside","outside"),    
                       limits=NULL,
                       weightClusterSize=1.0){
  
  n1<-dim(M)[1]
  if(is.list(clu)) {
    n<-sapply(clu, length)
  }  else{
    n<-length(clu)
  }
  if(sum(n)!=n1) stop("The length of clu and dimension of M does not match!")  
  diagonal<-match.arg(diagonal)
  limitType<-match.arg(limitType)  
  if(is.null(weights)){
    weights<-M
    weights[]<-1
  } else if(any(dim(weights)!=dim(M))) stop("Weights have wrong dim!")
  w<-weights
  if(is.null(uWeights)){
    uWeights<-rep(1.0, n1)
  }
  if(length(uWeights)!=n1) stop("uWeights has wrong length!")
  
  
  nMode<-ifelse(is.list(clu),length(clu),1)
  
  if(nMode>1){
    tmN<-sapply(clu,length)
    clu<-lapply(clu,function(x)as.integer(factor(x)))
    tmNclu<-sapply(clu,max)
    for(iMode in 2:nMode){
      clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
    }
    clu<-unlist(clu)
  } else{
    clu<-as.integer(factor(clu))
    tmNclu<-max(clu)
    tmN<-length(clu)
  }
  clu <- clu - 1
  if(length(dim(M))==2) M<-array(M,dim=c(dim(M),1))
  if(length(dim(w))==2) w<-array(w,dim=c(dim(w),1))
  
  if(is.null(limits)){
    bordersMatLower <- bordersMatUpper <- bordersSeperateLower <- bordersSeperateUpper<-NULL
    if(limitType!="none"){
      limitType<-"none"
      warning("limitType is set to 'none' as limits are NULL!")
    }      
  } else {
    if(diagonal %in% c("ignore","same")){
      bordersSeperateLower <- bordersSeperateUpper
      if(is.list(limits)){
        limits<-limits[[1]]
      }
      if(!is.array(limits)) stop("'limits' must be specified as an array!")
      dl<-dim(limits)
      if(length(dl)==3 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1:2],1,dl[3]))
      if(dim(limits)!=4) stop("'limits' has wrong dimmensions (see help for correct dimmensions)")
      
      if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
        stop("'limits' has wrong dimmensions (see help for correct dimmensions)")
      } else{
        bordersMatLower <- limits[,,,1] 
        bordersMatUpper <- limits[,,,2]
      }
    } else {
      if(is.list(limits) & length(limits)==2){
        diagLimits<-limits[[2]]
        limits<-limits[[1]]
      } else stop("If diagonal is 'seperate', limits must be a list of length 2")
      if(!is.array(limits)) stop("First element of 'limits' must be specified as an array!")
      dl<-dim(limits)
      if(length(dl)==3 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1:2],1,dl[3]))
      if(dim(limits)!=4) stop("First element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
        stop("First element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      } else{
        bordersMatLower <- limits[,,,1] 
        bordersMatUpper <- limits[,,,2]
      }
      
      if(!is.array(diagLimits)) stop("Second element of 'limits' must be specified as an array!")
      dl<-dim(diagLimits)
      if(length(dl)==2 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1],1,dl[2]))
      if(dim(diagLimits)!=3) stop("Second element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      if(all(dim(diagLimits)!=c(sum(tmNclu),dim(M)[3],2))){
        stop("Second element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      } else{
        bordersSeperateLower <- diagLimits[,,,1] 
        bordersSeperateUpper <- diagLimits[,,,2]
      }
    }
  }
  
  if(diagonal == "ignore")for(i in 1:dim(w)[3]){
    diag(w[,,i])<-0
  }
  w<-w*findEmptySuperbocks(M,n = n)
  w<-w/mean(w[w>0])
  uWeights<-uWeights/mean(uWeights[uWeights>0])
  weightClusterSize<-as.double(weightClusterSize)
  
  res<-.critFunction(M=M, clu=clu, weights=w, uWeights=uWeights, dimensions=sum(tmNclu), n=n, weightClusterSize=weightClusterSize, diagonal = diagonal, sBorders = limitType, bordersMatLower = bordersMatLower, bordersMatUpper = bordersMatUpper, bordersSeperateLower = bordersSeperateLower, bordersSeperateUpper = bordersSeperateUpper)
  
  # wByHB<-blockmodeling::funByBlocks(w,clu=rep(1:length(n),times=n),ignore.diag=FALSE, FUN=sum)
  # k<-length(unique(clu))
  # parByHB<-findActiveParam(M, n, k, na.rm=TRUE)
  # ICLpen<- sum(unclass(parByHB*log(wByHB)))+sum((k-1)*log(n))
  # ICL<- -res - 1/2*ICLpen

  return(ICL(M = M,clu = clu,weights = w,n = n,err=res))
  
  # res<-list(M=M, clu=clu, IM=IM, err=err, best=list(list(M=M, clu=clu, IM=IM)))
  # return(res)
}


# not exported
ICL<-function(M, clu, weights, n, err, ll){
  if(missing(err)) err<- -ll
  w<-weights
  wByHB<-blockmodeling::funByBlocks(w,clu=rep(1:length(n),times=n),ignore.diag=FALSE, FUN=sum)
  k<-length(unique(clu))
  parByHB<-findActiveParam(M, n, k, na.rm=TRUE)
  ICLpen<- sum(unclass(parByHB*log(wByHB)))+sum((k-1)*log(n))
  ICL<- -err - 1/2*ICLpen
  return(ICL)
}
