#' Relocation Heuristic for Generalized Structural Balance
#'
#' @description This function runs relocation heuristic for generalized structural balance on an \eqn{N x N} asymmetric matrix. The main diagonal is ignored.
#' @param A An \eqn{N x N} signed network matrix.
#' @param C The number of clusters (\eqn{1 < C < N}, where \eqn{N} is the number of nodes).
#' @param TLIMIT A desired time limit.
#' @return The function returns the following:
#' \itemize{
#' \item \code{obj} - the Doreian & Mrvar objective value;
#' \item \code{P} - \eqn{N}-dimensional vector of cluser assignements; and
#' \item \code{restarts} - the number of restarts within the time limit.
#' }
#' @examples
#' # Load the friendship network.
#' data("hafriend")
#'
#' # Run relocation heuristic for generalized structural balance.
#' res <- rhgsbt(A = hafriend, C = 3, TLIMIT = 1)
#'
#'# See the results.
#'res
#'
#' # Plot the network, in a matrix format, in line with the obtained partition P.
#' # The function plotMat is from blockmodeling package.
#' # plotMat(hafriend, clu = res$P)
#' @author Michael Brusco
#' @references
#' Brusco, M. J., & Doreian, P. (2019). Partitioning signed networks using relocation heuristics, tabu search, and variable neighborhood search. Social Networks, 56, 70-80. https://doi.org/10.1016/j.socnet.2018.08.007
#'
#' Doreian, P., & Mrvar, A. (2009). Partitioning signed social networks. Social Networks, 31, 1-11. http://dx.doi.org/10.1016/j.socnet.2008.08.001
#'
#' Doreian, P., & Mrvar, A. (1996). A partitioning approach to structural balance. Social Networks, 18, 149-168. https://doi.org/10.1016/0378-8733(95)00259-6

rhgsbt = function(A,C,TLIMIT) {
	N = dim(A)[1]
	OBJVAL = 0
  NREPS = 0
	EB <- matrix(0, nrow = N, ncol = 1)
	res <- .Fortran("rhgsbtf",as.integer(N),as.integer(C),as.double(TLIMIT),as.double(OBJVAL),as.integer(A),as.integer(EB),as.integer(NREPS))
	P <- res[[6]]
	#P <- matrix(p, nrow = N, ncol = 1)
	obj <- res[[4]]
	restarts <- res[[7]]
	return(list(P=P, obj=obj, restarts=restarts))
}
