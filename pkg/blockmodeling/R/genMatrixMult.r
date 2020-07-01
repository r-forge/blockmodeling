#' Generalized matrix multiplication
#' 
#' Computes a generalized matrix multiplication, where sum and product functions (elemet-wise and summary functions) can be replaced by arbitrary functions.
#' 
# #' @usage genMatrixMult(A, B, FUNelement = "*", FUNsummary = sum)
#' 
#' @param A The first matrix.
#' @param B The second matrix.
#' @param FUNelement Element-wise operator.
#' @param FUNsummary Summary function.
#'
#' @return A character vector or matrix.
#'
#' @examples
#' # Operations can be anything
#' x <- matrix(letters[1:8], ncol = 2)
#' y <- matrix(1:10, nrow = 2)
#' 
#' genMatrixMult(x, y, FUNelement = paste,
#' FUNsummary = function(x) paste(x, collapse = "|"))
#' 
#' # Binary logic
#' set.seed(1)
#' x <- matrix(rbinom(8, size = 1, prob = 0.5) == 1, ncol = 2)
#' y <- matrix(rbinom(10, size = 1, prob = 0.5) == 1, nrow = 2)
#' genMatrixMult(x, y, FUNelement = "*", FUNsummary = any)
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @seealso \code{\link{matmult}}
#' @keywords array algebra
#' 
#' @export

genMatrixMult<-function(A,B,FUNelement="*", FUNsummary=sum){
	if(dim(A)[2]!=dim(B)[1]) stop("incompatible dimmensions")
	n1<-dim(A)[1]
	n2<-dim(B)[2]
	X<-matrix(NA,nrow=n1,ncol=n2)
	dimnames(X)=list(dimnames(A)[[1]],dimnames(B)[[2]])
	for(i1 in 1:n1){
		for(i2 in 1:n2){
			X[i1,i2]<-FUNsummary(do.call(FUNelement,list(A[i1,],B[,i2])))
		}
	}
	return(X)
}
