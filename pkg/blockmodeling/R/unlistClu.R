#' @encoding UTF-8
#' @title Function for "unlisting" a partition. 
#' 
#' @description
#' Essentially, if the argument is a list (otherwise function just returns its argument), the function calls unlist on it. Before it, it however makes sure that names from different elements of the list to not repeat.  The opposite of \code{\link{splitClu}}. The \code{n} argument of the \code{\link{splitClu}} is returned as an attribute. If \code{renumber=TRUE} (default), it is practically identical to unlistCluInt.
#' 
#' @param clu A list representing a partition of units from different sets. Each element should be a partition for one set. 
#' @param renumber If \code{TRUE} (default), are renumbered so that they are 1:"total number of clusters". If any cluster "ID" is present in more than one set of units (one partition, one element of the list), this is done even if \code{renumber = FALSE}.
#'
#' 
#' @return A vector representing a partition. It also has an attribute \code{n} with the number of units that were in each set.
#' 
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' 
#' @seealso \code{\link{clu}}, \code{\link{splitClu}}, \code{\link{unlistCluInt}}
#'
#'
#' @examples
#' n <- c(8,8) 
#' cluList <- c(rep(1:2, times = c(3, 5)), rep(5:6, times = c(3, 5)))
#' unlistClu(clu = clu)
#' unlistClu(clu = clu, renumber = FALSE)
#'  
#' @keywords manip
#' @export


unlistClu<-function(clu, renumber=FALSE){
  if(!is.list(clu)){
	attr(clu, "n")<-length(clu)
	return(clu)
  }
  n<-sapply(clu, length)
  uniqueList<-lapply(clu, unique)
  k<-sapply(uniqueList,length)
  uniqueVec<-unlist(uniqueList)
  if(length(uniqueVec)!=length(unique(uniqueVec))) renumber<-TRUE
  if(renumber){
    clu<-lapply(clu, function(x)as.integer(factor(x)))
	for(i in 2:length(clu)){
      clu[[i]]<-clu[[i]]+sum(k[1:(i-1)])
	}
  }
  cluVec<-unlist(clu)
  attr(cluVec,"n")<-n
  return(cluVec)
}
