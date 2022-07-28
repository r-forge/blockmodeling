library(StochBlockTest)
library(blockmodeling)
#### Testing linked networks
set.seed(2021)
IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
clu<-rep(1:2, each=20)
n<-length(clu)
nClu<-length(unique(clu))
M1<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n)
M2<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n)
M12<-diag(n)
nn<-c(n,n)
k<-c(2,2)
Ml<-matrix(0, nrow=sum(nn),ncol=sum(nn))
Ml[1:n,1:n]<-M1
Ml[n+1:n,n+1:n]<-M2
Ml[n+1:n, 1:n]<-M12
plotMat(Ml)

ICLStochBlock(M=Ml,clu = list(clu,clu))

resMlbad<-stochBlock(M=Ml, clu = list(rep(1:2,n/2),rep(1:2,n/2)))
resMlOpt<-stochBlock(M=Ml, clu = list(clu,clu))
plot(resMlOpt)
crand(clu(resMlbad), clu(resMlOpt))

resMl<-stochBlockORP(M=Ml, k=k, n=nn, rep=100, nCores = 0)
plot(resMl)

resKMint<-stochBlockKMint(M=Ml, k=k, n=nn, nCores = 0)
plot(resKMint)
