#' @rdname Pajek
#' 
#' @description \code{savevector} - Saves a vector into a  Pajek ".clu" filename.
#'
#' @param v A vector.
#' 
#' @export

"savevector" <-
structure(function(v,filename){
if(length(grep(pattern="w32",x=version["os"]))){
	eol<-"\n"
}else{eol<-"\r\n"}
cat(paste(c(paste("*Vertices",length(v)), v),collapse=eol),file = filename)
}
, comment = "Save vector to file that can be read by Pajek")
