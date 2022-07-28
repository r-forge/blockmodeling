library(StochBlockTest)
library(blockmodeling)
packageVersion("StochBlockTest")
packageDate("StochBlockTest")


k<-2
blockSizes<-rep(20,k)
IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
if(any(dim(IM)!=c(k,k))) stop("invalid dimmensions")

set.seed(2021)
clu<-rep(1:k, times=blockSizes)
n<-length(clu)
M<-matrix(rbinom(n*n,1,IM[clu,clu]),ncol=n, nrow=n)
#clu<-clu-1
diag(M)<-0
plotMat(M)

cluK2<-rep(1:2,length.out=n)
resOptK2<-stochBlock(M,clu=cluK2)
resOptK2$ICL
resOptK2a1F<-stochBlock(M,clu=cluK2,addOne = FALSE)
resOptK2a1F$ICL
crand(clu(resOptK2),clu(resOptK2a1F))


cluK3<-rep(1:3,length.out=n)
resOptK3<-stochBlock(M,clu=cluK3)
resOptK3$ICL
#debugonce(stochBlock)
resOptK3a1F<-stochBlock(M,clu=cluK3,addOne = FALSE)
resOptK3a1F$ICL
crand(clu(resOptK3),clu(resOptK3a1F))
plot(resOptK3)
plot(resOptK3a1F)


resORP<-stochBlockORP(M,k=2, rep=100,nCores = 0,useParLapply = FALSE, return.all = TRUE)
resORP$ICL

resORPa1F<-stochBlockORP(M,k=2, rep=100,nCores = 0,useParLapply = FALSE, return.all = TRUE, addOne=FALSE)
resORPa1F$ICL
crand(clu(resORP),clu(resORPa1F))

resKMint<-stochBlockKMint(M,k=2,n=40, perm = 0,nCores = 0,useParLapply = FALSE, return.all = TRUE)
resKMint$ICL

resKMintList<-lapply(2:5, function(k)stochBlockKMint(M,k=k,n=40, perm = 0,nCores = 0,useParLapply = FALSE, return.all = TRUE))
sapply(resKMintList, function(x)x$ICL)


resKMintListPerm10<-lapply(2:5, function(k)stochBlockKMint(M,k=k,n=40, perm = 10,nCores = 0,useParLapply = FALSE, return.all = TRUE))
sapply(resKMintListPerm10, function(x)x$ICL)

resKMintListPerm100<-lapply(2:5, function(k)stochBlockKMint(M,k=k,n=40, perm = 100,nCores = 0,useParLapply = FALSE, return.all = TRUE))
sapply(resKMintListPerm100, function(x)x$ICL)




resORPList<-lapply(2:5, function(k)stochBlockORP(M,k=k,rep = 100,nCores = 0,useParLapply = FALSE, return.all = TRUE))


resORPa1FList<-lapply(2:5, function(k)stochBlockORP(M,k=k, rep = 100,nCores = 0,useParLapply = FALSE, return.all = TRUE, addOne=FALSE))


sapply(resORPList, function(x)x$ICL)
sapply(resORPa1FList, function(x)x$ICL)


vICL<-c(ICLStochBlock(M,clu=rep(1,nrow(M))),sapply(resORPList, function(x)x$ICL) )

vErr<-c(llStochBlock(M,rep(1,nrow(M))), sapply(resORPList, err))

plot(vICL,type="o")
plot(vErr,type="o")
-vICL-vErr

plot(resORP)
plot(resORPk3)
plot(resORPk4)


tmp<-cbind(ORP=sapply(resORPList, function(x)x$ICL),
      KMint=sapply(resKMintList, function(x)x$ICL),
      KMintPerm10=sapply(resKMintListPerm10, function(x)x$ICL),
      KMintPerm100=sapply(resKMintListPerm100, function(x)x$ICL))
plotMat(tmp-max(tmp))
plotMat(matrix(-100.001*(0:24), ncol=5))
debugonce(plotMat)
