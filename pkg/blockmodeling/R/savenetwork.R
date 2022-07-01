#' @rdname Pajek
#' 
#' @description \code{savenetwork} - Saves a matrix into a  Pajek ".net" filename.
#'
#' @param n A matrix representing the network.
#' @param twomode 1 for one-mode networks and 2 for two-mode networks. Default sets the argument to 1 for square matrices and to 2 for others.
#' @param symetric If \code{TRUE}, only the lower part of the matrix is used and the values are interpreted as "Edges", not "Arcs".
#' @import Matrix
#' @importFrom utils write.table
#' 
#' @export

"savenetwork" <-
structure(function(n,filename,twomode="default",symetric=NULL){
	if(length(grep(pattern="w32",x=version["os"]))){
		eol<-"\n"
	}else{eol<-"\r\n"}
	rowNames<-rownames(n)
	colNames<-colnames(n)
	if(requireNamespace("Matrix")){
		if(inherits(n,"mat")) n<-unclass(n)
		n <- as(n,"dgTMatrix")
		useMatrix<-TRUE
	}else{
		pack<-attr(class(n),"package")
		if(!(is.null(pack))&&pack=="Matrix") stop("The supplied object needs Matrix packege, but the package is not available.")
		useMatrix<-FALSE
	}
	if(dim(n)[1]!=dim(n)[2]){
		twomode<-2
	}else if(twomode=="default")twomode<-1
	if(is.null(symetric))if(twomode==1){
		if(useMatrix){symetric<-all(n==Matrix::t(n))
		}else symetric<-all(n==t(n))
	} else symetric<-FALSE
	pack<-attr("package",class(n))
	
	if ((dim(n)[1] == dim(n)[2]) & (twomode!=2)){
	  cat(paste("*Vertices",dim(n)[1]),eol, file = filename);
	  cat(paste(seq(1,length=dim(n)[1]),' "',rowNames,'"',eol,sep=""), file = filename,append=TRUE);
	  if(useMatrix){
	  	nDf<-as.data.frame(attributes(n)[c("i","j","x")])
	  	nDf[,c("i","j")]<-nDf[,c("i","j")]+1
	  	if(symetric){
			cat("*Edges",eol, file = filename,append=TRUE)	  		
			nDf<-nDf[nDf$i<=nDf$j,]
			write.table(nDf[,],eol=eol,file=filename,row.names = FALSE,col.names = FALSE,append=TRUE)
		} else {
			cat("*Arcs",eol, file = filename,append=TRUE)
			write.table(nDf[,],eol=eol,file=filename,row.names = FALSE,col.names = FALSE,append=TRUE)
		}
	  }else{
		  if(symetric){
			cat("*Edges",eol, file = filename,append=TRUE)
			   for (i in 1:dim(n)[1]) {
			     for (j in 1:(i)) {
			       if (n[i,j]!=0) {cat(paste(i,j,n[i,j],eol),file = filename,append=TRUE)}
			     }
			  } 
		  }else{
			   cat("*Arcs",eol, file = filename,append=TRUE);
			   for (i in 1:dim(n)[1]) {
			     for (j in 1:dim(n)[2]) {
			       if (n[i,j]!=0) {cat(paste(i,j,n[i,j],eol),file = filename,append=TRUE)}
			     }
			} 
		  }
	  } 
	}else { 
	  cat(paste("*Vertices",sum(dim(n)),dim(n)[1]),eol, file = filename);
	  cat(paste(1:dim(n)[1],' "',rowNames,'"',eol,sep=""), file = filename,append=TRUE);
	  cat(paste(seq(dim(n)[1]+1,length=dim(n)[2]),' "',colNames,'"',eol,sep=""), file = filename,append=TRUE);
	  cat("*Edges",eol, file = filename,append=TRUE);
	  if(useMatrix){
	  	nDf<-as.data.frame(attributes(n)[c("i","j","x")])
	  	nDf[,c("i","j")]<-nDf[,c("i","j")]+1
	  	nDf$j<-nDf$j+dim(n)[1]
	  	write.table(nDf[,],eol=eol,file=filename,row.names = FALSE,col.names = FALSE,append=TRUE)	  	
	  }else{
	   for (i in 1:dim(n)[1]) {
	     for (j in 1:dim(n)[2]) {
	       if (n[i,j]!=0) {cat(paste(i,j+dim(n)[1],n[i,j],eol),file = filename,append=TRUE)}
	     }
	   } 
	  }
	} 

}
, comment = "Save matrix to file that can be read by Pajek (as *Arcs)")
