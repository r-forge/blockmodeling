"loadvector" <-
structure(function(filename){
  i=0
  repeat {
    n<-read.table(file=filename,nrows=1,as.is=TRUE,skip = i)
    if(substr(n[1],1,1)!="%") break
    print(paste(n,collapse=" "))
    i=i+1
  }
  vv<-read.table(file=filename,skip=1+i,as.is=TRUE)
  if (dim(vv)[2]==1)
    vv<-vv[[1]]
  vv
}
, comment = "Load vector(s) from file that was produced by Pajek")
