#' @encoding UTF-8
#' @title Orders the partition so that mean values of \code{fun} applied to columns (if \code{funWay=2}, default), rows (if \code{funWay=1}) or both (if \code{funWay=c(1,2)}) is decreasing by clusters.
#' 
#' @description
#' Orders the partition so that mean values of \code{fun} applied to columns (if \code{funWay=2}, default), rows (if \code{funWay=1}) or both (if \code{funWay=c(1,2)}) is decreasing by clusters. The function can be used on the results of \code{\link{critFunC}}, \code{\link{optRandomParC}} or similar, or matrix and a partition can be supplied. It should also work on multirelational and lined networks.
#'
#' @param x A result of \code{\link{critFunC}}, \code{\link{optRandomParC}} or similar (something containing M (matrix) and clu (partition)) or a matrix (or array for multirelational networks).
#' @param clu A partition - a vector or a list of vectors/partitions. It must be supplied only if \code{x} is a matrix or array.
#' @param fun A function used to summarize rows or columns. \code{sum} by default.
#' @param funWay In which "way" should \code{fun} be appluied - to columns (if \code{funWay=2}, default), rows (if \code{funWay=1}) or both (if \code{funWay=c(1,2)})
#' @param nn The numbers of untis by sets of units. In principle, the function should determin this automatically.
#' @param returnList Logical. Should the partition be returned in form of a list (for lined networks only). \code{TRUE} by default.
#' @param scale Only used in case of multirelational networks. Should relations be scaled (\code{TRUE} by default) before summation. It can also be a vector of weights by relations.
#' @return An ordered partition. In an attribute ("reorder"). the information on how things were reordered.
#' @seealso \code{\link{clu}}
#' @export
orderClu<-function(x, clu=NULL,  fun=sum, funWay=2, nn=NULL, returnList=TRUE, scale=TRUE){
  if(inherits(x,c("check.these.par", "crit.fun", "critFun", "opt.more.par", "opt.more.par.mode", "opt.par", "opt.par.mode", "optMorePar", "optMoreParMode", "optPar", "optParMode"))){
    tclu<- clu(x)
    M<-x$M
    if(is.null(nn))nn<-x$initial.param$initial.param$n
  } else{
    M<-x
    if(is.null(clu)) stop("If x does not contain partition (clu), this must be supplied!")
    tclu<-clu
  }
  if(is.null(nn)&is.list(clu))nn<-sapply(clu,length)
  if(length(dim(M))>2){
    if(isFALSE(scale)){
      #do nothing
    }else if(isTRUE(scale)){
      myScale<-function(x)(x-mean(x))/stats::sd(x)
      for(i in 1:dim(M)[3])M[,,i]<-myScale(M[,,i])
    } else if(length(scale)==dim(M)[3]){
      for(i in 1:dim(M)[3])M[,,i]<-scale[i]*(M[,,i])
    }
    M<-apply(M,1:2, sum)
  }
  if(!is.null(nn)){
    rAll<-NULL
    m<-length(nn)
    tclu<- by(tclu, INDICES = rep(1:m, times=nn), FUN=c)
    k<-sapply(tclu,function(x)length(unique(x)))
    tcluAll<-NULL
    nCum<-cumsum(c(0,nn))
    kCum<-cumsum(c(0,k))
    for(i in 1:m){
      ids<-(nCum[i]+1):nCum[i+1]
      itclu<-tclu[[i]]
      iM<-M[ids, ids]
      crit<-unclass(by(data = apply(iM,funWay[1],fun, na.rm=TRUE),itclu,FUN = mean))
      if(length(funWay)==2) crit<-crit+unclass(by(data = apply(iM,funWay[2],fun, na.rm=TRUE),itclu,FUN = mean))
      r<-rank(-crit)+kCum[i]
      itclu<-r[as.character(itclu)]
      attr(itclu,"reorder")<-r
      rAll<-c(rAll,r)
      tcluAll<-c(tcluAll, list(itclu))
    } 
    if(!returnList) tcluAll<-unlist(tcluAll)
    attr(tcluAll,"reorder")<-rAll
    return(tcluAll)    
  }else{
    crit<-unclass(by(data = apply(M,funWay[1],fun, na.rm=TRUE),tclu,FUN = mean))
    if(length(funWay)==2) crit<-crit+unclass(by(data = apply(M,funWay[2],fun, na.rm=TRUE),tclu,FUN = mean))
    r<-rank(-crit)
    tclu<-r[as.character(tclu)]
    attr(tclu,"reorder")<-r
    return(tclu)
  }
}