library(StochBlockTest)

MlNet<-readRDS("alesPrimer.RDS")
n<-rep(nrow(MlNet)/4,4)
modes<-rep(1:4,times=n)
plotMat(MlNet, clu=modes)
k<-rep(3,4)

set.seed(2022)
resMlKMint<-stochBlockKMint(M=MlNet, k=k, n=n, nCores = 1, perm=9)
plot(resMlKMint)
