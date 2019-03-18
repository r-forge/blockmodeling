library(kmBlock)
library(blockmodelingTest)


k<-3
n<-20
IM<-outer(1:k,1:k)
clu<-sort(rep_len(1:k,length.out = n))
M<-IM[clu,clu]  
diag(M)<-100

nRep<-20

tmp<-system.time(test60sKmBlockC<-kmBlock::testORPtimeLim(M, k=k, blockFun = kmBlock::kmBlockC, timeLim = 60))

attr(test60sKmBlockC,which = "totalTime")<-tmp
test60sKmBlockC$nStarts




tmp<-system.time(test60sKmBlockR<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::kmBlock, timeLim = 60))

attr(test60sKmBlockR,which = "totalTime")<-tmp
test60sKmBlockR$nStarts



tmp<-system.time(test60sHomCmoveAndChange<-kmBlock::testORPtimeLim(M, k=k, blockFun = blockmodelingTest::optParC, timeLim = 0.000000000000000000000000001,approach="hom", blocks="com"))

attr(test60sHomCmoveAndChange,which = "totalTime")<-tmp
test60sHomCmoveAndChange$nStarts



homCmoveAndChangeTime<-system.time(res<-blockmodelingTest::optRandomParC(M,k=k,rep=nRep,approach="hom", blocks="com"))

homCjustMoveTime<-system.time(res<-blockmodelingTest::optRandomParC(M,k=k,rep=nRep,approach="hom", blocks="com",justChange=TRUE))

rbind(kmBlockORPCtime,kmBlockORPtime,homCmoveAndChangeTime,homCjustMoveTime)

#                       user.self sys.self elapsed user.child sys.child
# kmBlockORPCtime            0.90     0.00    0.91         NA        NA
# kmBlockORPtime             6.40     0.13    7.61         NA        NA
# homCmoveAndChangeTime    610.67     1.48  633.00         NA        NA
# homCjustMoveTime         389.79     0.19  392.47         NA        NA


