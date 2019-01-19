library(kmBlock)
library(blockmodelingTest)
M<-outer(rep(1:5, each=2),rep(1:5, each=2))
diag(M)<-100
M<-array(c(M,M*10), dim=c(10,10,2))
clu<-rep(0:4, each=2)
n<-10
nClu<-5

meanByBlocks(M, clu,nClu, diag="same",n )
meanByBlocks(M, clu,nClu, diag="ignore",n )
meanByBlocks(M, clu,nClu, diag="seperate",n )

fun.by.blocks(M[,,1],clu=clu)==meanByBlocks(M, clu,nClu, diag="ignore",n )$meansByBlocs[,,1]
fun.by.blocks(M[,,1],clu=clu,ignore.diag = FALSE)==meanByBlocks(M, clu,nClu, diag="same",n )$meansByBlocs[,,1]


#set weights to have the same dimensions as M and set al values to 1 (for testing)
W<-M
W[,,]<-1.0
M[5:10,5:10,]<-0
critFunction( M, clu, W, sum(nClu),n )
kmBlock::kmBlock(M,clu,W,sum(n),sum(nClu))

clu2<-c(0:4, 0:4)
critFunction( M, clu2, W, sum(nClu),n )
kmBlock::kmBlock(M,clu2,W,sum(n),sum(nClu))

clu3<-sample(c(0:4, 0:4))
critFunction( M, clu3, W, sum(nClu),n )
kmBlock::kmBlock(M,clu3,W,sum(n),sum(nClu))

n<-c(4,6)
nClu<-c(2,3)
clu4<-c(rep(0:1,times=2),rep(2:4,times=2))
critFunction( M[,,1,drop=FALSE], clu4, W, sum(nClu),n )
kmBlock::kmBlock(M[,,1,drop=FALSE],clu4,W,n,nClu)
blockmodelingTest::kmBlock(M[,,1],clu = clu4)[c("clu","err")]

set.seed(1)
clu4<-c(sample(rep(0:1,times=2)),sample(rep(2:4,times=2)))
critFunction( M[,,1,drop=FALSE], clu4, W, sum(nClu),n )
kmBlock::kmBlock(M[,,1,drop=FALSE],clu4,W,n,nClu)
blockmodelingTest::kmBlock(M[,,1],clu = clu4)[c("clu","err")]

n<-c(4,6)
nClu<-c(3,5)
clu5<-c(sample(c(0:2,2)),sample(c(3:7,7)))
critFunction( M[,,1,drop=FALSE], clu5, W, sum(nClu),n )
(tmp<-kmBlock::kmBlock(M[,,1,drop=FALSE],clu5,W,n,nClu))
blockmodelingTest::kmBlock(M[,,1],clu = clu5)[c("clu","err")]

critFunC(M, clu=tmp$bestClu,approaches = "hom",homFun = "ss", blocks=c("com"))$err
critFunction( M[,,1,drop=FALSE], tmp$bestClu, W, sum(nClu),n )
