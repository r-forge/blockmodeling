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

#set inf borders
bordersLower<-array(-Inf, dim=dim(M))
bordersUpper<-array(Inf, dim=dim(M))

bordersSeperateLower = matrix( -Inf, sum(nClu), dim(M)[3] )
bordersSeperateUpper = matrix( Inf, sum(nClu), dim(M)[3] )

#set weights to have the same dimensions as M and set al values to 1 (for testing)
W<-M
W[,,]<-1.0
M[5:10,5:10,]<-0
critFunction( M, clu, W, sum(nClu),n )
kmBlock::kmBlock(M,clu,W,sum(n),sum(nClu))

bordersLower<-array(-Inf, dim=c(sum(nClu),sum(nClu),dim(M)[3]))
bordersUpper<-array(Inf, dim=c(sum(nClu),sum(nClu),dim(M)[3]))

bordersSeperateLower = matrix( -Inf, sum(nClu), dim(M)[3])
bordersSeperateUpper = matrix( Inf, sum(nClu), dim(M)[3])

kmBlock::kmBlockC(M,clu,weights = W,  diagonal =  "seperate")

lim<- list(array(c(bordersLower,bordersUpper), dim=c(dim(bordersLower),2)), array(c(bordersSeperateLower,bordersSeperateUpper),dim=c(dim(bordersSeperateUpper),2)))
kmBlock::kmBlockC(M,clu,weights = W,  diagonal =  "seperate", limits =lim)

res<-kmBlock::kmBlockORPC(M,k=5,rep=100,weights = W,  diagonal =  "seperate")
system.time(res<-kmBlock::kmBlockORPC(M[,,1],k=5,rep=100,  diagonal =  "seperate"))

system.time(res<-blockmodelingTest::kmBlockORP(M[,,1],k=5,rep=100))

system.time(res<-blockmodelingTest::optRandomParC(M[,,1],k=5,rep=100,approach="hom", blocks="com"))


clu2<-c(0:4, 0:4)
critFunction( M, clu2, W, sum(nClu),n )
kmBlock::kmBlock(M,clu2,W,sum(n),sum(nClu))



clu3<-sample(c(0:4, 0:4))
critFunction( M, clu3, W, sum(nClu),n, "seperate", TRUE, bordersLower, bordersUpper, bordersSeperateLower, bordersSeperateUpper )
kmBlock::kmBlock(M,clu3,W,sum(n),sum(nClu), "seperate", TRUE, bordersLower, bordersUpper, bordersSeperateLower, bordersSeperateUpper )


lim2<-lim
lim2[[1]][3:5,3:5,,]<-0
lim2[[2]][3:5,,]<-0
clu3<-sample(c(0:4, 0:4))
clu3
kmBlock::kmBlockC(M,clu3,weights = W,  diagonal =  "seperate", limits =lim2)
kmBlock::kmBlockC(M,clu3,weights = W,  diagonal =  "seperate", limits =lim)
kmBlock::kmBlockC(M,clu3,weights = W,  diagonal =  "seperate")


lim3<-lim2
lim3[[1]][1:2,1:5,,1]<-rbind(1:5,(1:5)*2)
lim3[[1]][3:5,1:2,,1]<-cbind(3:5,(3:5)*2)
lim3[[1]][,,2,1]<-2*lim3[[1]][,,1,1]
lim3[[2]][1:2,,1]<-c(100,100,1000,1000)

for(i in 1:1000){
  clu3<-sample(c(0:4, 0:4))
  tmp1<-kmBlock::kmBlockC(M,clu3,weights = W,  diagonal =  "seperate", limits =lim3)
  tmp2<-kmBlock::kmBlockC(M,clu3,weights = W,  diagonal =  "seperate", limits =lim)
  tmp3<-kmBlock::kmBlockC(M,clu3,weights = W,  diagonal =  "seperate")
  
  if((tmp1$bestCf)==0) {
    print(tmp1)
    print(tmp2)
    print(tmp3)
    print(clu3)
    print(i)
    break
  }
}


n<-c(4,6)
nClu<-c(2,3)
clu4<-c(rep(0:1,times=2),rep(2:4,times=2))

bordersSeperateLower = matrix( -Inf, sum(nClu), dim(M)[3] )
bordersSeperateUpper = matrix( Inf, sum(nClu), dim(M)[3] )

critFunction( M[,,1,drop=FALSE], clu4, W, sum(nClu),n, "ignore", TRUE, bordersLower, bordersUpper )
kmBlock::kmBlock(M[,,1,drop=FALSE],clu4,W,n,nClu, "ignore", TRUE, bordersLower, bordersUpper )
blockmodelingTest::kmBlock(M[,,1],clu = clu4)[c("clu","err")]

set.seed(1)
clu4<-c(sample(rep(0:1,times=2)),sample(rep(2:4,times=2)))
critFunction( M[,,1,drop=FALSE], clu4, W, sum(nClu),n, "seperate", TRUE, bordersLower, bordersUpper, bordersSeperateLower, bordersSeperateUpper )
kmBlock::kmBlock(M[,,1,drop=FALSE],clu4,W,n,nClu, "seperate", TRUE, bordersLower, bordersUpper, bordersSeperateLower, bordersSeperateUpper )
blockmodelingTest::kmBlock(M[,,1],clu = clu4)[c("clu","err")]

n<-c(4,6)
nClu<-c(3,5)
clu5<-c(sample(c(0:2,2)),sample(c(3:7,7)))

bordersSeperateLower = matrix( -Inf, sum(nClu), dim(M)[3] )
bordersSeperateUpper = matrix( Inf, sum(nClu), dim(M)[3] )

critFunction( M[,,1,drop=FALSE], clu5, W, sum(nClu),n )
(tmp<-kmBlock::kmBlock(M[,,1,drop=FALSE],clu5,W,n,nClu))
blockmodelingTest::kmBlock(M[,,1],clu = clu5)[c("clu","err")]

critFunC(M, clu=tmp$bestClu,approaches = "hom",homFun = "ss", blocks=c("com"))$err
critFunction( M[,,1,drop=FALSE], tmp$bestClu, W, sum(nClu),n )
