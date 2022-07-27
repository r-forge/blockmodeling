#' @encoding UTF-8
#' @title Unlist a partition.
#' 
#' @description
#' It is used to convert a partition by sets into a single "simple" partition. Simple partition is a partition of only one set, that is a vector where units with the same value are considered to belong to the same cluster. The partitions by sets is a list, where each element of a list is a "simple" partition that corresponds to one set. The function first converts all elements of the lists to integers, that makes sure that each set uses different integers and on the end uses unlist function on such list.
#'
#' @param clu A partition by sets, that is a list of "simple" partitions.
#' @return The unlisted partition - one vector containing only integers.
#' @seealso \code{\link{clu}},  \code{\link{splitClu}}, \code{\link{unlistClu}}
#' @examples
#' cluList<-list(c("a","b","a"),c("b","c","b","c"))
#' unlistCluInt(cluList)
#' 
#' cluList<-list(c(1,1,1,2,2,2),c(1,1,1,2,2,2,3,3))
#' unlistCluInt(cluList)
#' @export
unlistCluInt<-function(clu){
  if(!is.list(clu)){
    warning("Clu must be a list! The orginal argument is returned")
  } else {
    clu<-lapply(clu,function(x)as.integer(as.factor(x)))
    nUnitsInRCclu<-lapply(clu,function(x)as.integer(table(x)))
    tmNclu<-sapply(clu,max)
    for(iMode in 2:length(clu)){
      clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
    }
    clu<-unlist(clu)
  }
  return(clu)
}