library(kmBlock)
library(blockmodelingTest)


k<-3
n<-20
IM<-outer(1:k,1:k)
clu<-sort(rep_len(1:k,length.out = n))
M<-IM[clu,clu]  
diag(M)<-100

nRep<-20

tmp<-system.time(test60sKmBlockC<-kmBlock::testORPtimeLim(M, k=k, blockFun = kmBlock::kmBlockC, timeLim = 5))

attr(test60sKmBlockC,which = "totalTime")<-tmp
test60sKmBlockC$nStarts




tmp<-system.time(test60sKmBlockR<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::kmBlock, timeLim = 5))

attr(test60sKmBlockR,which = "totalTime")<-tmp
test60sKmBlockR$nStarts



tmp<-system.time(test60sHomCmoveAndChange<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::optParC, timeLim = 5,approach="hom", blocks="com"))

attr(test60sHomCmoveAndChange,which = "totalTime")<-tmp
test60sHomCmoveAndChange$nStarts


tmp<-system.time(test60sHomCjustMove<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::optParC, timeLim = 5,approach="hom", blocks="com",justMove=TRUE))
attr(test60sHomCjustMove,which = "totalTime")<-tmp
test60sHomCjustMove$nStarts



tmp<-system.time(test60sHomCmoveAndChangeMulti<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::optParC, timeLim = 5,approach="hom", blocks="com", useMulti = TRUE))

attr(test60sHomCmoveAndChangeMulti,which = "totalTime")<-tmp
test60sHomCmoveAndChangeMulti$nStarts


tmp<-system.time(test60sHomCjustMoveMulti<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::optParC, timeLim = 5,approach="hom", blocks="com",justMove=TRUE, useMulti = TRUE))
attr(test60sHomCjustMoveMulti,which = "totalTime")<-tmp
test60sHomCjustMoveMulti$nStarts



seed<-sample(1000000,1)
set.seed(seed)
homCmoveAndChangeTime<-system.time(homCmoveAndChangeRes<-blockmodelingTest::optRandomParC(M,k=k,rep=nRep,approach="hom", blocks="com",useMulti=FALSE))
set.seed(seed)
homCjustMoveTime<-system.time(homCjustMoveRes<-blockmodelingTest::optRandomParC(M,k=k,rep=nRep,approach="hom", blocks="com",justChange=TRUE,useMulti=FALSE))
homCmoveAndChangeRes$err
homCjustMoveRes$err

rbind(kmBlockORPCtime,kmBlockORPtime,homCmoveAndChangeTime,homCjustMoveTime)

#                       user.self sys.self elapsed user.child sys.child
# kmBlockORPCtime            0.90     0.00    0.91         NA        NA
# kmBlockORPtime             6.40     0.13    7.61         NA        NA
# homCmoveAndChangeTime    610.67     1.48  633.00         NA        NA
# homCjustMoveTime         389.79     0.19  392.47         NA        NA


