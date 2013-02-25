"loadnetwork3" <-
function(filename,useSparseMatrix=NULL,minN=50){
  trim.trailing <- function (x) sub("\\s+$", "", x)
  rLines<-readLines(con=filename)
  nl<-length(rLines)
  ind.stars<-which(regexpr(pattern="*", text=rLines,fixed=TRUE)>0)
  nstars<-length(ind.stars)
  stars<-rLines[ind.stars]
  stars<-trim.trailing(stars)
  if(ind.stars[1]!=1){
    print(paste(rLines[1:(ind.stars[1]-1)],collapse="\n"))
  }
  rm(rLines)  
  n<-read.table(file=filename,skip=ind.stars[1]-1, nrows=1)
  print(paste(n,collapse=" "))
  if(length(n)==2){
    n<-as.numeric(n[2])
    vnames<-rep(as.character(NA),n)
    vnamesTab<-read.table(file=filename,skip=1,nrows=ind.stars[2]-ind.stars[1],as.is =TRUE)
    vnames[vnamesTab[,1]]<-vnamesTab[,2]
    if(all(is.na(vnames))){
        vnames<-NULL
    } else vnames[is.na(vnames)]<-""
    if(is.null(useSparseMatrix)){
        useSparseMatrix<- n>=50
    }
    
    if(useSparseMatrix){
    	if(require(Matrix)){
    		M<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    	}else{
    		warning("Matrix package is not installed. Ordanary (dense) matrices will be used instead of sparse onse")
    		M<-matrix(0,nrow=n,ncol=n)
    	}
    }else{
	M<-matrix(0,nrow=n,ncol=n)       
    }
    
    
    if(useSparseMatrix){
    	if(require(Matrix)){
    		M<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    	}else{
        	M<-matrix(0,nrow=n,ncol=n)
        	warning("Matrix package is not installed. Ordanary (dense) matrices will be used instead of sparse onse")
        }
    } else{
	M<-matrix(0,nrow=n,ncol=n)           
    }
    for(i in 2:nstars){
      nrows<-ifelse(i==nstars,-1,ind.stars[i+1]-ind.stars[i]-1)
      ties<-read.table(file=filename,skip=ind.stars[i],nrows=nrows)
      ncols<-dim(ties)[2]
      if(ncols==2){
      	ties<-cbind(ties,1)
      } else if(ncols>3){
      	ties<-ties[,1:3]
      }
      ties<-apply(ties,2,as.numeric)
      if(tolower(stars[i])=="*arcs"){
        M[ties[,1:2]]<-ties[,3]
      } else if(tolower(stars[i])=="*edges"){
        M[ties[,1:2]]<-ties[,3]
        M[ties[,2:1]]<-ties[,3]
      }
    }
    dimnames(M)<-list(vnames,vnames)
  } else{
      n12<-as.numeric(n[2])
      n1<-as.numeric(n[3])
      n2<-n12-n1
      vnames1<-read.table(file=filename,skip=1,nrows=n12)[,2]
    vnames<-read.table(file=filename,skip=1,nrows=n12,as.is =TRUE)[,2]
    if(all(is.na(vnames))){
        vnames<-NULL
    } else vnames[is.na(vnames)]<-""
    rLines<-readLines(con=filename)
    nl<-length(rLines)
    ind.stars<-which(regexpr(pattern="*", text=rLines,fixed=TRUE)>0)
    nstars<-length(ind.stars)
    stars<-rLines[ind.stars]
    rm(rLines)
    if(is.null(useSparseMatrix)){
        useSparseMatrix<- n12>50
    }    
    if(useSparseMatrix){
    	if(require(Matrix)){
    		M<-Matrix(0,nrow=n12,ncol=n12,sparse=TRUE)
    	}else{
    		warning("Matrix package is not installed. Ordanary (dense) matrices will be used instead of sparse onse")
    		M<-matrix(0,nrow=n12,ncol=n12)
    	}
    }else{
	M<-matrix(0,nrow=n12,ncol=n12)       
    }
    for(i in 2:nstars){
      nrows<-ifelse(i==nstars,-1,ind.stars[i+1]-ind.stars[i]-1)
      ties<-read.table(file=filename,skip=ind.stars[i],nrows=nrows)
      ncols<-dim(ties)[2]
      if(ncols==2){
      	ties<-cbind(ties,1)
      } else if(ncols>3){
      	ties<-ties[,1:3]
      }
      ties<-apply(ties,2,as.numeric)
      M[ties[,1:2]]<-ties[,3]
      M[ties[,2:1]]<-ties[,3]
    }
    dimnames(M)<-list(vnames,vnames)
  M<-M[1:n1,(n1+1):n12]    
  }
  return(M)
}

