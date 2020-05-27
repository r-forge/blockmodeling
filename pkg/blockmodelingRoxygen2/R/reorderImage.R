#' Reordering an image matrix of the blockmodel (or an error matrix based on new and old partition
#' 
#' Reorders an image matrix of the blockmodel (or an error matrix based on new and old partition.
#' The partitions should be the same, except that classes can have different labels.
#' It is useful when we want to have a different order of classes in figures and then also in image matrices. Currently it is only suitable for one-mode blockmodels.
#'
# #' @usage reorderImage(IM, oldClu, newClu)
#' 
#' @param IM An image or error matrix.
#' @param oldClu Old partition.
#' @param newClu New partition, the same as the old one except for class labeles.
#'
#' @return Reorder matrix (rows and columns are reordred).
#' 
#' @references
#' \enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#'
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#'
#' @author Ales Ziberna
#' @seealso \code{\link{critFunC}}, \code{\link{plot.mat}}, \code{\link{clu}}, \code{\link{IM}}, \code{\link{err}}
#' @keywords manip
#' 
#' @export

reorderImage<-function(IM,oldClu,newClu){
if(crand2(oldClu,newClu)!=1)stop("Old and new clu's are not compatibale (crand index is not 1)!\n")
newOrder<-which(table(oldClu,newClu)>0,arr.ind=TRUE)[,1]
return(IM[newOrder,newOrder])
}
