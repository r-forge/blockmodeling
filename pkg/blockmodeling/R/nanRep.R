#' @encoding UTF-8
#' @title Replaces NaN values by the speficied values (0 by default)
#'
#' @param x A vector or similar where the NaNs are to be replaced.
#' @param rep A value that should replace the NaNs (0 by default).
#' @return x with NaNs replaced by rep.
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' @keywords manip
#' @export
nanRep<-function(x, rep=0){
  x[is.nan(x)]<-rep
  return(x)
}