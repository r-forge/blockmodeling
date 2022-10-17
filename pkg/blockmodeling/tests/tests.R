library(blockmodeling)
nCores<-1
clu <-  c(1, 2, 1, 2, 1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 4, 3, 4, 2, 2, 3)
data(baker)

set.seed(2022)
res<-optRandomParC(baker>0,k=4, approaches = "bin", blocks = c("nul","com"),rep = 100, nCores = nCores)
plot(res)
print(res)


resSS<-optRandomParC(baker>0,k=4, approaches = "hom", blocks = c("nul","com"),rep = 100, nCores = nCores)
plot(resSS)
print(resSS)

resBll<-optRandomParC(baker>0,k=4, approaches = "hom", blocks = c("nul","com"),rep = 100, nCores = nCores, homFun="bll")
plot(resBll)
print(resBll)

if(requireNamespace("StochBlockTest")){
  StochBlockTest::llStochBlock(baker>0, clu=clu(resBll), addOne = FALSE, diagonal = "seperate")
  resSB<-StochBlockTest::stochBlockORP(baker>0,k = 4, rep = 100, addOne = FALSE, diagonal = "seperate")
  err(resSB)
  plot(resSB)
  crand(clu(resSB),clu(resBll))
}


tmp<-critFunC(baker>0, clu=clu, approaches = "hom", blocks = c("nul","rre"),homFun="bll", mulReg = TRUE)
plot(tmp)
tmp[["IM"]][1,,]
tmp[["EM"]][1,,]


tmp<-critFunC(baker>0, clu=clu, approaches = "hom", blocks = c("nul","com"),homFun="bll", mulReg = TRUE)
tmp$EM[1,,]
tmp$err
plot(tmp)
critFunC(baker>0, clu=clu, approaches = "hom", blocks = c("nul","com"),homFun="bll", mulReg = TRUE, diag=2)$err

#if(requireNamespace("StochBlockTest")) StochBlockTest::llStochBlock(baker>0, clu=clu,addOne = FALSE)

clu2L<-list(rep(1:2, each=5),rep(1:2, each=5))
tmp<-critFunC(baker>0, clu=clu2L, approaches = "hom", blocks = c("nul","com"),homFun="bll", mulReg = TRUE)

tmp$EM[1,,]
tmp$err

