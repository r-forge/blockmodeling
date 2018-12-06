library(kmBlock)
M<-outer(rep(1:5, each=2),rep(1:5, each=2))
diag(M)<-100
M<-array(c(M,M*10), dim=c(10,10,2))
clu<-rep(0:4, each=2)
n<-10
nClu<-5

meanByBlocks(M, clu,nClu, diag="same")
meanByBlocks(M, clu,nClu, diag="ignore")
meanByBlocks(M, clu,nClu, diag="seperate")


#criterialFunction( M, clu, weights, sum(nClu) )