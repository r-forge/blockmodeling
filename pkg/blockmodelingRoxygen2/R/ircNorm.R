#' Function for iterated row and column normalization of valued matrices
#' 
#' The aim is to obtain a matrix with row and column sums equal to 1.
#' This is achieved by iterating row and column normalization. This is usually not possible if any row or column has only 1 non-zero cell.
#' 
# #' @usage ircNorm(M, eps = 10^-12, maxiter = 1000)
#'
#' @param M A non-negative valued matrix to be normalized.
#' @param eps The maximum allows squared deviation of a row or column's maximum from 1 (if not exactly 0).
#' Also, if the all deviations in two consequtive iterations are smaller, the process is terminated.
#' @param maxiter Maximum number of iterations. If reached, the process is terminated and the current solution returned.
#'
#' @return Normalized matrix.
#'
#' @examples
#' A <- matrix(runif(100), ncol = 10)
#' A # A non-normalized matrix with different row and column sums.
#' apply(A, 1, sum)
#' apply(A, 2, sum)
#' A.norm <- ircNorm(A)
#' A.norm # Normalized matrix with all row and column sums approximately 1.
#' apply(A.norm, 1, sum)
#' apply(A.norm, 2, sum)
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @keywords manip

ircNorm<-function(M,eps=10^-12,maxiter=1000){
  diffM<-function(M){
    max(c(1-apply(M,1,sum)[apply(M,1,sum)>0],1-apply(M,2,sum)[apply(M,2,sum)>0])^2)
  }
  side<-1	
  i=0
  tmpM<-list(M,M)
  while(diffM(M)>eps){
    i=i+1
    sums<-apply(M, side, sum)
    sums[sums==0]<-1
    M<-sweep(M, side, sums,FUN="/")
    if(max(c(M-tmpM[[side]])^2)<eps) {
      warning("The covergence (in terms of row/column sums beeing equal to 1 not possible. Covergence reach in terms of stability of matrix after each transformation.")
      break
    }
    tmpM[[side]]<-M
    side<-3-side	
    if(i>=maxiter){
      warning("Maximum number of itrerations (",maxiter,") reached, convergence not achieved.\n")
      break
    }
  }	
  M<-(tmpM[[1]]+tmpM[[2]])/2
  return(M)
}