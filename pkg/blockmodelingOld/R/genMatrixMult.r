genMatrixMult<-function(A,B,FUNelement="*", FUNsummary=sum){
	if(dim(A)[2]!=dim(B)[1]) stop("incompatible dimmensions")
	n1<-dim(A)[1]
	n2<-dim(B)[2]
	X<-matrix(NA,nrow=n1,ncol=n2)
	dimnames(X)=list(dimnames(A)[[1]],dimnames(B)[[2]])
	for(i1 in 1:n1){
		for(i2 in 1:n2){
			X[i1,i2]<-FUNsummary(do.call(FUNelement,list(A[i1,],B[,i2])))
		}
	}
	return(X)
}
