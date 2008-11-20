"savenetwork" <-
structure(function(n,filename,twomode="default",symetric=NULL){
	rowNames<-rownames(n)
	colNames<-colnames(n)
	if(require(Matrix)){
		if(class(n)=="mat") n<-unclass(n)
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
	  write(paste("*Vertices",dim(n)[1]), file = filename);
	  write(paste(seq(1,length=dim(n)[1]),' "',rowNames,'"',sep=""), file = filename,append=TRUE);
	  if(useMatrix){
	  	nDf<-as.data.frame(attributes(n)[c("i","j","x")])
	  	nDf[,c("i","j")]<-nDf[,c("i","j")]+1
	  	if(symetric){
			write("*Edges", file = filename,append=TRUE)	  		
			nDf<-nDf[nDf$i<=nDf$j,]
			write.table(nDf[,],file=filename,row.names = FALSE,col.names = FALSE,append=TRUE)
		} else {
			write("*Arcs", file = filename,append=TRUE)
			write.table(nDf[,],file=filename,row.names = FALSE,col.names = FALSE,append=TRUE)
		}
	  }else{
		  if(symetric){
			write("*Edges", file = filename,append=TRUE)
			   for (i in 1:dim(n)[1]) {
			     for (j in 1:(i)) {
			       if (n[i,j]!=0) {write(paste(i,j,n[i,j]),file = filename,append=TRUE)}
			     }
			  } 
		  }else{
			   write("*Arcs", file = filename,append=TRUE);
			   for (i in 1:dim(n)[1]) {
			     for (j in 1:dim(n)[2]) {
			       if (n[i,j]!=0) {write(paste(i,j,n[i,j]),file = filename,append=TRUE)}
			     }
			} 
		  }
	  } 
	}else { 
	  write(paste("*Vertices",sum(dim(n)),dim(n)[1]), file = filename);
	  write(paste(1:dim(n)[1],' "',rowNames,'"',sep=""), file = filename,append=TRUE);
	  write(paste(seq(dim(n)[1]+1,length=dim(n)[2]),' "',colNames,'"',sep=""), file = filename,append=TRUE);
	  write("*Edges", file = filename,append=TRUE);
	  if(useMatrix){
	  	nDf<-as.data.frame(attributes(n)[c("i","j","x")])
	  	nDf[,c("i","j")]<-nDf[,c("i","j")]+1
	  	nDf$j<-nDf$j+dim(n)[1]
	  	write.table(nDf[,],file=filename,row.names = FALSE,col.names = FALSE,append=TRUE)	  	
	  }else{
	   for (i in 1:dim(n)[1]) {
	     for (j in 1:dim(n)[2]) {
	       if (n[i,j]!=0) {write(paste(i,j+dim(n)[1],n[i,j]),file = filename,append=TRUE)}
	     }
	   } 
	  }
	} 

}
, comment = "Save matrix to file that can be read by Pajek (as *Arcs)")
