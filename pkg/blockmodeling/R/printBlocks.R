#' @encoding UTF-8
#' @title Nice printing of the \code{blocks} parameter as used in \code{\link{optRandomParC}} and \code{\link{critFunC}}.
#'  
#' @param blocks \code{blocks} parameter as used in \code{\link{optRandomParC}} and \code{\link{critFunC}}.
#' 
#'
#' @return 
#' Used for side effects (printing)
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' @seealso \code{\link{optRandomParC}}, \code{\link{critFunC}}
#' @keywords print
#' 
#' @export
printBlocks<-function(blocks){
  B<-blocks
  if(is.vector(B)){
    if(is.list(B)){
      for(i in 1:length(B)){
        cat("Relation ",i,":", sep="")
        cat(B[[i]])
      }
    } else cat(B,"\n")
  } else{
    if(length(dim(B))==2){
      print(data.frame(B,check.names = FALSE))
    } else if(length(dim(B))==3){
      print(data.frame(apply(B,2:3,function(x)paste(na.omit(x),collapse=",")),check.names = FALSE))
    } else if(length(dim(B))==4){
	  if(dim(B)[2]==1){
		printBlocks(B[,1,,])
	  } else for(i in 1:dim(B)[2]){
        cat("Relation",i,"\n")
        print(data.frame(apply(B[,i,,],2:3,function(x)paste(na.omit(x),collapse=",")),check.names = FALSE))
      }
    }
  }
}