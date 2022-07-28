library(StochBlockTest)
library(blockmodeling)
#source("C:/Users/zibernaa/OneDrive - Univerza v Ljubljani/Multilevel Blockmodeling/stochBlockMay2021.R")
#source("C:/Users/zibernaa/OneDrive - Univerza v Ljubljani/Multilevel Blockmodeling/stochBlockMar2021.R")

set.seed(2021)
IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
clu<-rep(1:2, each=5)
n<-length(clu)
nClu<-length(unique(clu))
M<-matrix(rbinom(100,1,IM[clu,clu]),ncol=n, nrow=n)
clu<-clu-1
diag(M)<-0
plotMat(M)
M<-array(c(M,1-M), dim=c(n,n,2))
M[,,1]




W<-M
W[,,]<-1.0
W[1:5,6:10,]<-0
W1<-M
W1[,,]<-1.0

StochBlockTest:::.critFunction(M[,,1,drop=FALSE], clu = clu, weights = W1[,,1,drop=FALSE], dimensions = sum(nClu), n = n,weightClusterSize = 1,uWeights = rep(1,10))

StochBlockTest:::.critFunction(M[,,1,drop=FALSE], clu = clu, weights = W1[,,1,drop=FALSE], dimensions = sum(nClu), n = n,weightClusterSize = 1,uWeights = rep(1,10),addOne = FALSE)


#llStochBlockR(M=M[,,1],clu = clu,weights = W1[,,1])$err
StochBlockTest::llStochBlock(M=M[,,1],clu = clu)
ICLStochBlock(M=M[,,1],clu = clu)

ICLStochBlock(M=M[,,1],clu = clu,addOne = FALSE)
#debugonce(ICL)

ICLStochBlock(M=M,clu = clu)



#debugonce(StochBlockTest:::ICL)
ICLStochBlock(M=M[,,1],clu = rep(c(0,1),c(n-1,1)))


clu2<-clu
clu2[2]<-1
set.seed(2021)
#tmpR<-stochBlock(M = M[,,1],clu = clu2)
set.seed(2021)
tmpC<-StochBlockTest::stochBlock(M = M[,,1],clu = clu2)
tmpC$err
tmpC$ICL
#tmpR$err
clu(tmpC)
#clu(tmpR)

tmpCMr<-StochBlockTest::stochBlock(M = M,clu = clu2)
#debugonce(ICL)



#StochBlockTest:::.critFunction(M[,,1,drop=FALSE], clu = clu, weights = W1[,,1,drop=FALSE], dimensions = sum(nClu), n = n,weightClusterSize = 1,uWeights = 1)
#llStochBlockR(M=M[,,1],clu = clu,weights = W1[,,1])$err
StochBlockTest::llStochBlock(M=M[,,1],clu = clu)



