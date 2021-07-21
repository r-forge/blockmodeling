#' @encoding UTF-8
#' @title Create canonical partition and find unique canonical partitions in a list of partitions.
#' 
#' @description
#' It is used to convert any partition to a canonical partition. A canonical partition is a partition where the first unit is in cluster 1, the next unit that is not in cluster 1 in in cluster 2 and so on. So if we would take first appearances of clusters in the order they appear in the partition vector, we would get integers from 1 to the number of clusters.
#'
#' @param clu A partition - a vector or a list of vectors/partitions.
#' @param cluList A list of partitions(vectors).
#' @return For function \code{canClu} - a canonical partition or a list of such partitions.
#' For function \code{canCluUniqe} - A list of unique canonical partitions.
#' @seealso \code{\link{clu}}
#' @examples
#' clu<-c(3,2,2,3,1,2)
#' canClu(clu)
#' @export
canClu<-function(clu){
  if(!is.list(clu)){
    return(as.numeric(factor(clu,levels=unique(clu))))
  } else {
    lapply(clu, canClu)
  }
}


#' @rdname canClu
#' 
#' @export
canCluUniqe<-function(cluList){
  if(!is.list(cluList)){
    stop("cluList must be a list of partitions!")
  } else {
    uniqueClu<-NULL
    uniqueCluStr<-NULL
    cluList<-lapply(cluList, canClu)
    cluListStr<-sapply(cluList, paste, collapse=",")
    for(i in 1:length(cluList)){
      if(!(cluListStr[i]%in%uniqueCluStr)){
        uniqueClu<-c(uniqueClu,cluList[i])
        uniqueCluStr<-c(uniqueCluStr,cluListStr[i])
      }
    }
    return(uniqueClu)
  }
}