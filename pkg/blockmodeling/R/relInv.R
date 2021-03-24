#' @encoding UTF-8
#' @title Functions for computing "relative inverse" (\code{x[1]/x}).
#' 
#' @description
#' For a vector x, it computes \code{x[1]/x}. For \code{relInv2}, if certain elements of the result are not finite (e.g. if certain elements of x are 0), these elements are replaced with 0s.
#'
#' @param x A numeric vector. For \code{relInv} it should not contain 0s (while for  \code{relInv2} it can).
#' @return A vector computed as \code{x[1]/x}. For \code{relInv2}, if the non-finite elements are replaced with 0s.
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' @keywords manip
#' @export
relInv<-function(x)x[1]/x


#' @rdname relInv
#' 
#' @export
relInv2<-function(x){
  x<-x[1]/x
  x[!is.finite(x)]<-0
  x
}

