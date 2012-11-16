"loadnetwork2" <-
function(filename,useSparseMatrix=NULL,minN=50,safe=TRUE,closeFile=TRUE){
  if(is.character(filename)){
  	file<-file(description=filename,open="r")
  } else file<-filename
  while(TRUE){
	line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
	if(substr(line[1],start=1,stop=1)=="%") {print(paste(line,collapse=" "));next}
	n<-line
	break
  }
  if(length(n)==2){
    n<-as.numeric(n[2])
    if(safe){
		vnames<-rep(as.character(NA),n)
		while(TRUE){
			line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
			if(length(line)==0||sum(grep(patt="^ *$",x=as.character(line))==1)) break
			if(substr(line[1],start=1,stop=1)=="%") {print(paste(line,collapse=" "));next}
			if(substr(line[1],start=1,stop=1)=="*"){
				type=line[1]
				break
			}
			vnames[as.integer(line[1])]<-line[2]
    	}
    }else{
    	vnames<-read.table(file=file,nrows=n,as.is =TRUE)[,2]
    	type=""   	
    }
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
	        M<-matrix(0,nrow=n,ncol=n)
	        warning("Matrix package is not installed. Ordanary (dense) matrices will be used instead of sparse onse")    	
    	}
    }else{
        M<-matrix(0,nrow=n,ncol=n)
    }
    
    if(type=="*Matrix"){
        tmp<-read.table(file=file,nrows=n)
        tmp<-as.matrix(tmp)
        M[1:n,1:n]<-M
    } else while(TRUE){
    	line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
    	if(length(line)==0||sum(grep(patt="^ *$",x=as.character(line))==1)) break
    	if(substr(line[1],start=1,stop=1)=="%") {print(paste(line,collapse=" "));next}
    	if(substr(line[1],start=1,stop=1)=="*"){
    		type=line[1]
    		next
    	}else line<-as.double(line)
    	
    	if(type=="*Arcs"){
    		M[line[1],line[2]]<-line[3]
    	}else if(type=="*Edges") {
    		M[line[1],line[2]]<-line[3]
    		M[line[2],line[1]]<-line[3]
    	}	
    }
    dimnames(M)<-list(vnames,vnames)
  } else{
    n12<-as.numeric(n[2])
    n1<-as.numeric(n[3])
    n2<-n12-n1
    if(safe){
		vnames<-rep(as.character(NA),n12)
		while(TRUE){
			line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
			if(length(line)==0||sum(grep(patt="^ *$",x=as.character(line))==1)) break
			if(substr(line[1],start=1,stop=1)=="%") {print(paste(line,collapse=" "));next}
			if(substr(line[1],start=1,stop=1)=="*"){
				type=line[1]
				break
			}
			vnames[as.integer(line[1])]<-line[2]
    	}
    }else{
    	vnames<-read.table(file=file,nrows=n,as.is =TRUE)[,2]
    	type=""   	
    }

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
    } else {
	M<-matrix(0,nrow=n12,ncol=n12)       
    }
    if(type=="*Matrix"){
        tmp<-read.table(file=file,nrows=n1)
        tmp<-as.matrix(tmp)
        M[1:n1,(n1+1):n12]<-tmp
    } else while(TRUE){
    	line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
    	if(length(line)==0||sum(grep(patt="^ *$",x=as.character(line))==1)) break
    	if(substr(line[1],start=1,stop=1)=="%") {print(paste(line,collapse=" "));next}
    	if(substr(line[1],start=1,stop=1)=="*"){
    		type=line[1]
    		next
    	}else line<-as.double(line)
    	
	M[line[1],line[2]]<-line[3]
	M[line[2],line[1]]<-line[3]
    }
    dimnames(M)<-list(vnames,vnames)
    M<-M[1:n1,(n1+1):n12]    
  }
  if(closeFile) close(file)
  return(M)
}
