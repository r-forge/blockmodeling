"loadnetwork" <-
function(filename,useSparseMatrix=NULL,minN=50){
  trim.trailing <- function (x) sub("\\s+$", "", x)

  n<-c("%")
  i=0
  repeat {
    n<-read.table(file=filename,nrows=1,as.is=TRUE,skip = i)
    if(substr(n[1],1,1)!="%") break
    print(paste(n,collapse=" "))
    i=i+1
  }
  if(length(n)==2){
    n<-as.numeric(n[2])
    n
    rLines<-readLines(con=filename)
    nl<-length(rLines)
    #ind.stars<-which(regexpr(pattern="*", text=rLines,fixed=TRUE)>0)
    ind.stars<-which(substr(rLines,1,1)=="*")
    nstars<-length(ind.stars)
    stars<-rLines[ind.stars]
    stars<-trim.trailing(stars)
    rm(rLines)
    
    vnames1<-read.table(file=filename,skip=ind.stars[1],nrows=ind.stars[2]-ind.stars[1]-1,as.is =TRUE)
    vnames<-character(n)
    vnames[vnames1[,1]]<-vnames1[,2]
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
      if(stars[i]=="*Arcs"|stars[i]=="*arcs"){
        M[ties[,1:2]]<-ties[,3]
      } else if(stars[i]=="*Edges"|stars[i]=="*edges"){
        M[ties[,1:2]]<-ties[,3]
        M[ties[,2:1]]<-ties[,3]
      }
    }
    dimnames(M)<-list(vnames,vnames)
  } else{
      n12<-as.numeric(n[2])
      n1<-as.numeric(n[3])
      n2<-n12-n1

    rLines<-readLines(con=filename)
    nl<-length(rLines)
    #ind.stars<-which(regexpr(pattern="*", text=rLines,fixed=TRUE)>0)
    ind.stars<-which(substr(rLines,1,1)=="*")
    nstars<-length(ind.stars)
    stars<-rLines[ind.stars]
    rm(rLines)
    
    vnames1<-read.table(file=filename,skip=ind.stars[1],nrows=ind.stars[2]-ind.stars[1]-1,as.is =TRUE)
    vnames<-character(n12)
    vnames[vnames1[,1]]<-vnames1[,2]
    if(all(is.na(vnames))){
      vnames<-NULL
    } else vnames[is.na(vnames)]<-""
    
    
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

