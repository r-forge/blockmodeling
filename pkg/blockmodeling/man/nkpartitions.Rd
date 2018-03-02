\name{nkpartitions}
\alias{nkpartitions}
\alias{nkpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for listing all possible partitions or just cunting the number of them.}
\description{
The function \code{nkpartitions} lists all possible partitions of n objects in to k clusters. The function \code{nkpar} only gives the number of such partitions.
}
\usage{
nkpartitions(n, k, exact = TRUE, print = FALSE)
nkpar(n, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of units/objects}
  \item{k}{Number of clusters/groups}
  \item{exact}{Search for partitions with exactly k or at most k clusters}
  \item{print}{print results as they are found?}
}
\value{
The matrix or number of possible partitions.
}
\author{Chris Andrews}
\examples{
n<-8 #if larger, the number of partitions increases dramaticaly,
     #as does if we increase the number of clusters
net<-matrix(NA,ncol=n,nrow=n)
clu<-rep(1:2,times=c(3,5))
tclu<-table(clu)
net[clu==1,clu==1]<-rnorm(n=tclu[1]*tclu[1],mean=0,sd=1)
net[clu==1,clu==2]<-rnorm(n=tclu[1]*tclu[2],mean=4,sd=1)
net[clu==2,clu==1]<-rnorm(n=tclu[2]*tclu[1],mean=0,sd=1)
net[clu==2,clu==2]<-rnorm(n=tclu[2]*tclu[2],mean=0,sd=1)

#computation of criterion function with the correct partition
nkpar(n=n, k=length(tclu)) #computing the number of partitions
all.par<-nkpartitions(n=n, k=length(tclu)) #forming the partitions
all.par<-lapply(apply(all.par,1,list),function(x)x[[1]])
# to make a list out of the matrix
res<-check.these.par(M=net,partitions=all.par,approach="ss",blocks="com")
plot(res) #we get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS

