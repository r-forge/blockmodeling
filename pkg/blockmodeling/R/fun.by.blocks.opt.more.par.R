#' @rdname funByBlocks
#' @param orderClu Should the partition be ordered before computing. \code{FALSE} by default. If \code{TRUE}, \code{\link{orderClu}} is used (using default arguments) to order the clusters in a partition in "decearsing" (see \code{\link{orderClu}} for interpretation) order.  If \code{TRUE}, \code{sortNames} is set to \code{FALSE}.
#' @export
"funByBlocks.optMorePar" <-
function(
  x,	#an object of class "optMorePar"
  which=1,	#which best solution/partition should be used
  orderClu=FALSE, #should the partition be ordered.
  sortNames =  NULL, 
  ...	#aditional parameters to function "funByBlocks"
){
	if(which>length(x$best)){
		which<-1
		warning("Only",length(x$best),"solutions exists. The first solution will be used.")
	}
	tclu<-clu(x,which=which)
	if(orderClu) {
		tclu<-orderClu(x=x$M, clu=tclu)
		if(is.null(sortNames))sortNames<-FALSE
	} else if(is.null(sortNames))sortNames<-TRUE
	
	funByBlocks(M=x$M, clu=tclu,...)
}

#' @rdname funByBlocks
#' @method funByBlocks opt.more.par
#' @export
funByBlocks.opt.more.par<-funByBlocks.optMorePar
