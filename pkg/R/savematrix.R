"savematrix" <-
structure(function(n,filename,twomode=1){
if ((dim(n)[1] == dim(n)[2]) & (twomode!=2))
{ 
  verNames<-rownames(n)
  verNamesTable<-table(verNames)
  if(max(verNamesTable)>1){
  	duplicateName<-names(which(verNamesTable>1))
  	for(i in duplicateName){
  		verNames[verNames==i]<-paste(i,1:verNamesTable[i],sep="")
  	}
  }
  write(paste("*Vertices",dim(n)[1]), file = filename);
  write(paste(seq(1,length=dim(n)[1]),' "',verNames,'"',sep=""), file = filename,append=TRUE);
  write("*Matrix", file = filename,append=TRUE);
  write(t(n),file = filename,ncolumns=dim(n)[1],append=TRUE) }
else
{ 
  verRowNames<-rownames(n)
  verRowNamesTable<-table(verRowNames)
  if(max(verRowNamesTable)>1){
  	duplicateRowName<-names(which(verRowNamesTable>1))
  	for(i in duplicateRowName){
  		verRowNames[verRowNames==i]<-paste(i,1:verRowNamesTable[i],sep="")
  	}
  }
  verColNames<-colnames(n)
  verColNamesTable<-table(verColNames)
  if(max(verColNamesTable)>1){
  	duplicateColName<-names(which(verColNamesTable>1))
  	for(i in duplicateColName){
  		verColNames[verColNames==i]<-paste(i,1:verColNamesTable[i],sep="")
  	}
  }
  write(paste("*Vertices",sum(dim(n)),dim(n)[1]), file = filename);
  write(paste(1:dim(n)[1],' "',verRowNames,'"',sep=""), file = filename,append=TRUE);
  write(paste(seq(dim(n)[1]+1,length=dim(n)[2]),' "',verColNames,'"',sep=""), file = filename,append=TRUE);
  write("*Matrix", file = filename, append=TRUE);
  write(t(n),file = filename, ncolumns=dim(n)[2],append=TRUE)} }
, comment = "Save matrix to file that can be read by Pajek (as *Matrix)")
