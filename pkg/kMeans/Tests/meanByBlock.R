library(kmBlock)
library(blockmodelingTest)
M<-outer(rep(1:5, each=2),rep(1:5, each=2))
diag(M)<-100
M<-array(c(M,M*10), dim=c(10,10,2))
clu<-rep(0:4, each=2)
n<-10
nClu<-5

meanByBlocks(M, clu,nClu, diag="same")
meanByBlocks(M, clu,nClu, diag="ignore")
meanByBlocks(M, clu,nClu, diag="seperate")

fun.by.blocks(M[,,1],clu=clu)==meanByBlocks(M, clu,nClu, diag="ignore")$meansByBlocs[,,1]
fun.by.blocks(M[,,1],clu=clu,ignore.diag = FALSE)==meanByBlocks(M, clu,nClu, diag="same")$meansByBlocs[,,1]


#set weights to have the same dimensions as M and set al values to 1 (for testing)
W<-M
W[,,]<-1.0

critFunction( M, clu, W, sum(nClu) )
