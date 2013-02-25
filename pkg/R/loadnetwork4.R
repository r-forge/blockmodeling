"loadnetwork4" <-
function(filename,useSparseMatrix=NULL,minN=50,fill=FALSE){
  sc<-scan(filename,what="raw",sep="\n")
  sc<-gsub(patt="\\",rep="/",x=sc,fixed=TRUE)
  first<-sapply(sc,substr,start=1,stop=1)
  sc<-sc[first!="%"]
  first<-first[first!="%"]
  stars<-which(first=="*")
  stars<-c(stars,"*end"=length(sc)+1)
  n<-as.numeric(strsplit(sc[1]," +")[[1]][-1])
  if(is.null(useSparseMatrix)){
    useSparseMatrix<- n[1]>=minN
  }
  if(length(n)==1){
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
    
    vnames<-rep(as.character(""),n)
    
    for(i in seq_along(stars)){
      #i<-1
      type<-strsplit(x=names(stars)[i],split=" +")[[1]][1]
      if(tolower(type)=="*vertices"){
        #vnames<-rep(as.character(NA),n)
        verNames<-sc[(stars[i]+1):(stars[i+1]-1)]
        verNames<-paste(verNames,collapse="\n")
        verNames<-read.table(text=verNames,as.is=TRUE,fill=fill)
        vnames[verNames[,1]]<-verNames[,2]
      } else if(tolower(type)%in%c("*arcs","*edges")){
        ties<-sc[(stars[i]+1):(stars[i+1]-1)]
        ties<-paste(ties,collapse="\n")
        ties<-read.table(text=ties)
        ncols<-dim(ties)[2]
        if(ncols==2){
          ties<-cbind(ties,1)
        } else if(ncols>3){
          ties<-ties[,1:3]
        }
        ties<-apply(ties,2,as.numeric)
        if(tolower(type)=="*arcs"){
          M[ties[,1:2]]<-ties[,3]
        } else if(tolower(type)=="*edges"){
          M[ties[,1:2]]<-ties[,3]
          M[ties[,2:1]]<-ties[,3]
        }
      }
      dimnames(M)<-list(vnames,vnames)      
      
    }
    
  } else{
    n12<-n[1]
    n1<-n[2]
    n2<-n12-n1
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
    
    
    vnames<-rep(as.character(""),n12)
    
    for(i in seq_along(stars)){
      #i<-1
      type<-strsplit(x=names(stars)[i],split=" +")[[1]][1]
      if(tolower(type)=="*vertices"){
        #vnames<-rep(as.character(NA),n12)
        verNames<-sc[(stars[i]+1):(stars[i+1]-1)]
        verNames<-paste(verNames,collapse="\n")
        verNames<-read.table(text=verNames,as.is=TRUE,fill=fill)
        vnames[verNames[,1]]<-verNames[,2]
      } else if(tolower(type)%in%c("*arcs","*edges")){
        ties<-sc[(stars[i]+1):(stars[i+1]-1)]
        ties<-paste(ties,collapse="\n")
        ties<-read.table(text=ties)
        ncols<-dim(ties)[2]
        if(ncols==2){
          ties<-cbind(ties,1)
        } else if(ncols>3){
          ties<-ties[,1:3]
        }
        ties<-apply(ties,2,as.numeric)
        if(tolower(type)=="*arcs"){
          M[ties[,1:2]]<-ties[,3]
        } else if(tolower(type)=="*edges"){
          M[ties[,1:2]]<-ties[,3]
          M[ties[,2:1]]<-ties[,3]
        }
      }
      dimnames(M)<-list(vnames,vnames)      
      
    }
    
    
    M<-M[1:n1,(n1+1):n12]    
  }

  return(M)
}