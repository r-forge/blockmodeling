library(kmBlock)
library(blockmodelingTest)


k<-5
n<-200
IM<-outer(1:k,1:k)
clu<-sort(rep_len(1:k,length.out = n))
M<-IM[clu,clu]  
diag(M)<-100

nRep<-20

kmBlockORPCtime<-system.time(res<-kmBlock::kmBlockORPC(M,k=k,rep=nRep,  diagonal =  "seperate"))

kmBlockORPtime<-system.time(res<-blockmodelingTest::kmBlockORP(M,k=k,rep=nRep))

homCmoveAndChangeTime<-system.time(res<-blockmodelingTest::optRandomParC(M,k=k,rep=nRep,approach="hom", blocks="com"))

homCjustMoveTime<-system.time(res<-blockmodelingTest::optRandomParC(M,k=k,rep=nRep,approach="hom", blocks="com",justChange=TRUE))

rbind(kmBlockORPCtime,kmBlockORPtime,homCmoveAndChangeTime,homCjustMoveTime)

#                       user.self sys.self elapsed user.child sys.child
# kmBlockORPCtime            0.90     0.00    0.91         NA        NA
# kmBlockORPtime             6.40     0.13    7.61         NA        NA
# homCmoveAndChangeTime    610.67     1.48  633.00         NA        NA
# homCjustMoveTime         389.79     0.19  392.47         NA        NA


