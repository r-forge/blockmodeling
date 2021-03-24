#' @encoding UTF-8
#' @title Functions creating a list of partitions based on a single partition and information on the number of units in each set.
#' 
#' @description
#' Function \code{splitClu} creates a list of partitions based on a single partition (\code{clu}) and information on the number of units in each set (\code{n}).
#' 
#' Function \code{splitCluRes} does the same but extracts the information from the result of (old versions of) functions \code{\link{critFunC}}, \code{\link{optParC}}, \code{\link{optRandomParC}} or similar (newer versions should already return a list of partitions in case they are used on networks with more sets of units.
#' 
#' @param clu A vector representing a partition of units from different sets. Result of some legacy code for \code{\link{optRandomParC}} or \code{\link{optParC}} or similar functions.
#' @param n A vector with number of units per set. The assuption is that the first \code{n[1]} elements of \code{clu} are for the first set, the second \code{n[2]} elements of \code{clu} are for the second set and so on. \code{sum(n)} must be equal to \code{length(clu)}.
#' @param res Result of  (old versions of) functions \code{\link{critFunC}}, \code{\link{optParC}}, \code{\link{optRandomParC}} or similar.
#' @param renumber If \code{TRUE} (default), dimensions that have only one level are dropped
#' (\code{drop} function is applied to the final result).
#' @param \dots Not used.
#'
#' 
#' @return A list of partitions if \code{clu}, one for each set of units. A single vector if only one set of units is present.
#' 
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' 
#' @seealso \code{\link{clu}}
#'
#'
#' @examples
#' n <- c(8,8) 
#' clu <- c(rep(1:2, times = c(3, 5)), rep(3:4, times = c(3, 5)))
#' splitClu(clu = clu,n = n)
#' splitClu(clu = clu,n = n,renumber = TRUE)
#'  
#' @keywords manip
#' @export


splitClu<-function(clu, n, renumber=FALSE){
  if(length(n)==1) return(clu)
  clu<-split(clu, rep(1:length(n),times=n))
  if(renumber){
    return(lapply(clu, function(x)as.integer(factor(x))))
  }else return(clu)
}


#' @rdname splitClu
#' 
#' @export
splitCluRes<-function(res, renumber=FALSE){
  clu<-clu(res)
  n<-res$initial.param$n
  if(is.list(clu)) return(clu)
  splitClu(clu,n,renumber=renumber)
}
