ircNorm<-function(M,eps=10^-12,maxiter=1000){
	diffM<-function(M){
		max(c(1-apply(M,1,sum)[apply(M,1,sum)>0],1-apply(M,2,sum)[apply(M,2,sum)>0])^2)
	}
	side<-1	
	i=0
	tmpM<-list(M,M)
	while(diffM(M)>eps){
		i=i+1
		sums<-apply(M, side, sum)
		sums[sums==0]<-1
		M<-sweep(M, side, sums,FUN="/")
		if(max(c(M-tmpM[[side]])^2)<eps) {
		    warning("The covergence (in terms of row/column sums beeing equal to 1 not possible. Covergence reach in terms of stability of matrix after each transformation.")
		    break
		}
		tmpM[[side]]<-M
		side<-3-side	
		if(i>=maxiter){
			warning("Maximum number of itrerations (",maxiter,") reached, convergence not achieved.\n")
			break
		}
	}	
	M<-(tmpM[[1]]+tmpM[[2]])/2
	return(M)
}