#' Two-Mode Blockmodeling (Structural Equivalence) Heuristic
#'
#' @description This function runs two-mode KL-medians for an \eqn{RO x CO} network matrix.
#' @param A An \eqn{RO x CO} two-mode network matrix.
#' @param RC The number of clusters for row objects (\eqn{1 < RC < RO}).
#' @param CC The number of clusters for column objects (\eqn{1 < CC < CO}).
#' @param TLIMIT A desired time limit.
#' @return The function returns the following:
#' \itemize{
#' \item \code{objval} - total number of inconsistencies;
#' \item \code{RP} - an \eqn{RO}-dimensional vector of row cluser assignements;
#' \item \code{RC} - an \eqn{RC}-dimensional vector of column cluser assignements;
#' \item \code{restarts} - the number of restarts within the time limit.
#' }
#' @examples
#' # Load the Galaskiewicz’ (1985) CEO/Club network data.
#' data("galask")
#'
#' # Run the two-mode blockmodeling heuristic procedure.
#' res <- tmklmed(galask, RC = 5, CC = 5, TLIMIT = 1)
#'
#' # See the results.
#' res
#'
#' # Plot the network, in a matrix format, in line with the obtained partitions.
#' # The function plotMat is from blockmodeling package.
#' # plotMat(galask, clu = list(res$RP, res$CP))
#' @author Michael Brusco
#' @references
#' Brusco, M. J., & Doreian, P. (2019). Partitioning signed networks using relocation heuristics, tabu search, and variable neighborhood search. Social Networks, 56, 70-80. https://doi.org/10.1016/j.socnet.2018.08.007
#'
#' Doreian, P., & Mrvar, A. (2009). Partitioning signed social networks. Social Networks, 31, 1-11. http://dx.doi.org/10.1016/j.socnet.2008.08.001
#'
#' Doreian, P., & Mrvar, A. (1996). A partitioning approach to structural balance. Social Networks, 18, 149-168. https://doi.org/10.1016/0378-8733(95)00259-6

tmklmed = function(A,RC,CC,TLIMIT) {
	RO = dim(A)[1]
	CO = dim(A)[2]
	GBEST = 0
      NREPS = 0
	GR <- matrix(0, nrow = RO, ncol = 1)
	GC <- matrix(0, nrow = CO, ncol = 1)
	res =.Fortran("tmklmedf",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.integer(A),as.integer(GR),as.integer(GC),as.integer(GBEST),as.integer(NREPS))
	RP <- unlist(res[[7]])
	CP <- res[[8]]
	objval <- res[[9]]
	restarts <- res[[10]]
	return(list(RP=RP, CP=CP, objval=objval, restarts=restarts))
}
