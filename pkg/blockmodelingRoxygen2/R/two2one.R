#' Two-mode network conversions
#' 
#' Converting  two mode networks from two to one mode matrix representation and vice versa.
#' If a two-mode matrix is converted into a one-mode matrix, the original two-mode matrix lies in the upper right corner of the one-mode matrix.
#' 
# #' @usage
# #' two2one(M, clu = NULL)
# #' one2two(M, clu = NULL)
#'
#' @param M A matrix representing the (usually valued) network.
#' @param clu A partition. Each unique value represents one cluster. This should be a list of two vectors, one for each mode.
#'
#' @return
#' Function returns list with the elements:
#' a two mode matrix of a the two mode network in its upper left corner.
#' \item{M}{The matrix.}
#' \item{clu}{The partition, in form appropriate for the mode of the matrix.}
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @seealso \code{\link{optParC}}, \code{\link{optParC}}, \code{\link{optRandomParC}}, \code{\link{plot.mat}}
#'
#' @examples
#' # Generating a simple network corresponding to the simple Sum of squares
#' # Structural equivalence with blockmodel:
#' # null com
#' # null null
#' n <- c(7, 13)
#' net <- matrix(NA, nrow = n[1], ncol = n[2])
#' clu <- list(rep(1:2, times = c(3, 4)), rep(1:2, times = c(5, 8)))
#' tclu <- lapply(clu, table)
#' net[clu[[1]] == 1, clu[[2]] == 1] <- rnorm(n = tclu[[1]][1] * tclu[[2]][1],
#'    mean = 0, sd = 1)
#' net[clu[[1]] == 1, clu[[2]] == 2] <- rnorm(n = tclu[[1]][1] * tclu[[2]][2],
#'    mean = 4, sd = 1)
#' net[clu[[1]] == 2, clu[[2]] == 1] <- rnorm(n = tclu[[1]][2] * tclu[[2]][1],
#'    mean = 4, sd = 1)
#' net[clu[[1]] == 2, clu[[2]] == 2] <- rnorm(n = tclu[[1]][2] * tclu[[2]][2],
#'    mean = 0, sd = 1)
#' plot.mat(net, clu = clu) # Two mode matrix of a two mode network
#'
#' # Converting to one mode network
#' M1 <- two2one(net)$M
#' plot.mat(M1, clu = two2one(net)$clu) # Plotting one mode matrix
#' # Converting one to two mode matrix and plotting
#' plot.mat(one2two(M1, clu = clu)$M, clu = clu)
#' 
#' @keywords cluster graphs
#' 
#' @export

"two2one" <-
function(M,clu=NULL){
	n1<-dim(M)[1]
	n2<-dim(M)[2]
	n<-n1+n2
	M1<-matrix(0,nrow=n,ncol=n)
	M1[1:n1,(n1+1):n]<-M
	dimnames(M1)<-list(unlist(dimnames(M)),unlist(dimnames(M)))
  if(!is.null(clu)) {
    clu<-lapply(clu,function(x)as.numeric(as.factor(x)))
    clu[[2]]<-clu[[2]]+max(clu[[1]])
    clu<-unlist(clu)
  }
 return(list(M=M1,clu=clu))
}

