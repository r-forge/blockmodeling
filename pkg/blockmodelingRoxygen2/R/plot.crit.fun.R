#' @rdname plotMat
#' @export
"plot.critFun" <-
function(
	x,#an "critFun" class object
	main=NULL,
	...	#aditional parameters to "plot.mat"
){
	if(is.null(main)) main <- deparse(substitute(x))
	plot.mat(x$M,clu=x$clu,IM=x$IM,main=main,...)
}


#' @rdname plotMat
#' @method plot crit.fun
#' @export
plot.crit.fun<-plot.critFun