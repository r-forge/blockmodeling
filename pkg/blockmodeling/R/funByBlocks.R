#' Computation of function values by blocks
#' 
#' Computes a value of a function over blocks of a matrix, defined by a partition.
#' 
#' @param x An object of suitable class or a matrix/array representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation.
#' The network can have one or more modes (different kinds of units with no ties among themselves.
#' If the network is not two-mode, the matrix must be square.
#' @param clu A partition. Each unique value represents one cluster.
#' If the network is one-mode, then this should be a vector, else a list of vectors, one for each mode.
#' @param M A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation.
#' The network can have one or more modes (different  kinds of units with no ties among themselves.
#' If the network is not two-mode, the matrix must be square.
#' @param ignore.diag Should the diagonal be ignored.
#' @param sortNames Should the rows and columns of the matrix be sorted based on their names.
#' @param FUN The function to be computed over the blocks.
#' @param which Which (if several) of the "best" solutions should be used.
#' @param \dots Further arguments to \code{funByBlocks.default}. 
#'
#' @return A numerical matrix of \code{FUN} values by blocks, induced by a partition \code{clu}.
#' 
#' @references \enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#'
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @seealso \code{\link{optRandomParC}}, \code{\link{optParC}}
#' @examples
#' n <- 8 # If larger, the number of partitions increases dramatically,
#' # as does if we increase the number of clusters
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(3, 5))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#' # Optimizing 10 random partitions with optRandomParC
#' res <- optRandomParC(M = net, k = 2, rep = 10, approaches = "hom", homFun = "ss", blocks = "com")
#' plot(res) # Hopefully we get the original partition
#' funByBlocks(res)
#' # Computing mean by blocks, ignoring the diagonal (default)
#' 
#' @keywords cluster math
#' 
#' @export
funByBlocks <-
function(x, ...) UseMethod("funByBlocks")

#' @rdname funByBlocks
#' @export
fun.by.blocks<-funByBlocks