\name{fun.by.blocks}
\alias{fun.by.blocks}
\alias{fun.by.blocks.default}
\alias{fun.by.blocks.opt.more.par}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of function values by blocks}
\description{
Computes a value of a functions over blocks of a matrix, defined by a partition.
}
\usage{
fun.by.blocks(x, ...)

## default S3 method:
fun.by.blocks(M, clu, ignore.daig = identical(ss(diag(M)) && !is.list(clu), 0), FUN = "mean")

## S3 method for class 'opt.more.par':
fun.by.blocks(x, which = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{clu}{A partition. Each unique value represents one cluster. If the nework is  one-mode, than this should be a vector, else a list of vectors, one for each mode}
  \item{ignore.daig}{Should the diagonal be ingored. }
  \item{FUN}{Function to be computed over the blocks}
  
  \item{x}{An object of an appropriate class}
  \item{which}{Which (if several) of the "best" solutions should be used}
  \item{\dots}{Further arguments to \code{fun.by.blocks.default} }

}

\value{
 A numerical matrix of \code{FUN} values by blocks, induced by a partition \code{clu}
}
\references{
ŽIBERNA, Aleš (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

ŽIBERNA, Aleš. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.
}

\author{Aleš Žiberna}
\seealso{\code{\link{opt.random.par}},\code{\link{opt.these.par}}}

\examples{
n<-8 #if larger, the number of partitions increases dramaticaly, as does if we increase the number of clusters
net<-matrix(NA,ncol=n,nrow=n)
clu<-rep(1:2,times=c(3,5))
tclu<-table(clu)
net[clu==1,clu==1]<-rnorm(n=tclu[1]*tclu[1],mean=0,sd=1)
net[clu==1,clu==2]<-rnorm(n=tclu[1]*tclu[2],mean=4,sd=1)
net[clu==2,clu==1]<-rnorm(n=tclu[2]*tclu[1],mean=0,sd=1)
net[clu==2,clu==2]<-rnorm(n=tclu[2]*tclu[2],mean=0,sd=1)

#we select a random parition and then optimise it

all.par<-nkpartitions(n=n, k=length(tclu)) #forming the partitions
all.par<-lapply(apply(all.par,1,list),function(x)x[[1]]) # to make a list out of the matrix

#optimizing 10 random chosen partitions which are suplied to opt.these.par
res<-opt.these.par(M=net,partitions=all.par[sample(1:length(all.par),size=10)],approach="ss",blocks="com")
plot(res) #Hopefully we get the original partition
fun.by.blocks(res) #computing mean by blocks, ignoring the diagonal (default)
res$best[[1]]$BM #the same result computed by opt.these.par when approach="ss" and blocks="com"
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
\keyword{ math }% at least one, from doc/KEYWORDS
