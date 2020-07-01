#' @encoding UTF-8
#' @title Recode
#' 
#' @description
#' Recodes values in a vector.
#' 
# #' @usage recode(x, oldcode = sort(unique(x)), newcode)
#' 
#' @param x A vector.
#' @param oldcode A vector of old codes.
#' @param newcode A vector of new codes.
#'
#' @return A recoded vector.
#'
#' @examples
#' x <- rep(1:3, times = 1:3)
#' newx <- recode(x, oldcode = 1:3, newcode = c("a", "b", "c"))
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' @keywords manip
#' 
#' @export

"recode" <-
function(x,oldcode=sort(unique(x)),newcode){
  if(length(oldcode)!=length(newcode))stop("The number of old and new codes do not match")
  newx<-x
  for(i in 1:length(oldcode)){
    newx[x==oldcode[i]]<-newcode[i]
  }
  return(newx)
}

