#' @encoding UTF-8
#' @title Comparing partitions on one or multiple sets of units
#' 
#' @description
#' Rand Index and Rand Index corrected/adjusted for chance for comparing partitions (Hubert & Arabie, 1985). The functions also support computing these indices on partitions on multiple sets (where a "combined" partition is a list of multiple partitions).
#' The names of the clusters do not matter.
#' 
#' 
# #' @usage
# #' rand(clu1, clu2, tab)
# #' crand(clu1, clu2, tab)
#' 
#' @param clu1 The first of the two partitions to be compared, given in the form of vectors, where for each unit a cluster membership is given. Alternatively, this can be a contingency table obtained as a \code{table(clu1, clu2)}. If a partition, \code{clu2} must also be provided. In case of multiple sets, this should be pa list of partitions.
#' @param clu2 If \code{clu1} is partition or a list of partitions, this must be a comaptible the second partition or list of partitions.
#' @param tab A contingency table obtained as a \code{table(clu1, clu2)}. This is included for back-compatibility reasons. If this is present, all other arguments are ignored.
#' @param multiSets How should we compute the index in case of multiple sets of unis (if \code{clu1} and  \code{clu2} are lists of partitions)? Possible values are "unlist" and "weight" (the default).
#' @param weights Weights to be used if \code{multiSets} is "weight". It can be "equal", "size" (default) or a numeric (non-negative) vector of the same length as the number of sets (the number of partitions in the list of partitions).
#' @param returnIndividual If \code{multiSets} is "weight", should the indices for individual sets be also returned.  If \code{TRUE}, the function returns a list instead of a single value. If the values is \code{"attr"} (the default), the indices by sets are given as an attribute \code{"bySets"} 
#'
#' @return The value of Rand Index (corrected/adjusted for chance) unless \code{multiSets="weight"} and \code{returnIndividual=FALSE}. In this case, a list with two items is return. The "global" index is in \code{global}, while the the indices by sets are in \code{bySets}.
#'
#' @references Hubert, L., & Arabie, P. (1985). Comparing Partitions. Journal of Classification, 2(1), 193-218.
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @keywords cluster
#' 
#' @export 

"rand" <-
function (clu1, clu2, tab) #Hubert & Arabie
{
  if(missing(tab))if(is.table(clu1)){
      tab<-clu1
    } else tab<-table(clu1,clu2)
    
    n <- sum(tab)
	1 + (sum(tab^2) - (sum(rowSums(tab)^2) + sum(colSums(tab)^2))/2)/choose(n, 2)
}


#' @rdname rand
#' 
#' @export

"crand" <-
  function (clu1,clu2, tab, multiSets=c("weights","unlist"), weights = c("size","equal"), returnIndividual="attr") #Hubert & Arabie
  {
    if(missing(tab)) if(is.table(clu1)){
      tab<-clu1
    } else {
      if(is.list(clu1)){
        if(!is.list(clu2)|(length(clu1)!=length(clu2))) stop("If clu1 is a list, clu2 must be a list of equal size!")
        multiSets<-match.arg(multiSets)
        if(multiSets=="unlist"){
          tab<-table(unlistCluInt(clu1),unlistCluInt(clu2))
        }else if(multiSets=="weights"){
          if(is.numeric(weights)){
            if(length(weights)!=length(clu1)) stop("Weigts must equal the number of sets")
          } else {
            weights<-match.arg(weights)
            if(weights=="equal") {
              weights<-rep(1,length(clu1))
            }  else if(weights=="size"){
              weights<-sapply(clu1,length)
            } else stop("Unexpected 'weights' argument!")
          }
          bySets<-mapply(crand2, clu1, clu2)
          global<-sum(bySets*weights)/sum(weights)
          if(returnIndividual=="attr"){
            attr(global,"bySets")<-bySets
            return(global)
          } else if(returnIndividual==TRUE){
            return(list(global=global, bySets=bySets)) 
          } else if(returnIndividual==FALSE){
            return(global)
          } else stop("Unexpected 'returnIndividual' argument!")
          
        } else stop("Unexpected 'multiSets' argument!")
        
      } else tab<-table(clu1,clu2)
    }

    n <- sum(tab)
    sum.ni2 <- sum(choose(rowSums(tab), 2))
    sum.nj2 <- sum(choose(colSums(tab), 2))
    E<- sum.ni2 * sum.nj2 / choose(n, 2)
    return((sum(choose(tab, 2)) - E)/((sum.ni2 + sum.nj2)/2 - E))
}



#' @rdname rand
#' 
#' @export
"rand2" <-
  function (clu1, clu2) #Hubert & Arabie
  {
#    .Deprecated("rand")
    tab<-table(clu1,clu2)
    
    n <- sum(tab)
    1 + (sum(tab^2) - (sum(rowSums(tab)^2) + sum(colSums(tab)^2))/2)/choose(n, 2)
  }



#' @rdname rand
#' 
#' @export

"crand2" <-
  function (clu1,clu2) #Hubert & Arabie
  {
#    .Deprecated("crand")
    tab<-table(clu1,clu2)
    
    n <- sum(tab)
    sum.ni2 <- sum(choose(rowSums(tab), 2))
    sum.nj2 <- sum(choose(colSums(tab), 2))
    E<- sum.ni2 * sum.nj2 / choose(n, 2)
    return((sum(choose(tab, 2)) - E)/((sum.ni2 + sum.nj2)/2 - E))
}

