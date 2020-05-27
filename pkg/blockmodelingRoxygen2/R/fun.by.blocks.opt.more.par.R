#' @rdname funByBlocks
#' @export
"funByBlocks.optMorePar" <-
function(
  x,	#an object of class "optMorePar"
  which=1,	#which best solution/partition should be used
  ...	#aditional parameters to function "funByBlocks"
){
	if(which>length(x$best)){
		which<-1
		warning("Only",length(x$best),"solutions exists. The first solution will be used.")
	}
	funByBlocks(M=x$M, clu=clu(x,which=which),...)
}

