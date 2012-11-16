\name{two2one}
\alias{two2one}
\alias{one2two}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Two-mode network conversions}
\description{
  Coverting two mode networks from two to one mode matrix representation and vice versa. If a two-mode matrix is converted in-to a one-mode matrix, the original two-mode matrix lies in the upper right corner of the one-mode matrix.
}
\usage{
two2one(M, clu = NULL)
one2two(M, clu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network.}
  \item{clu}{A partition. Each unique value represents one cluster. This should be a list of two vectors, one for each mode.}
}

\value{
Functions returns list with elemets:
a mode mode matrix with the two mode network in its upper left corner.
  \item{M}{The matrix}
  \item{clu}{The partition, in form appropriate for the mode of the matrix}

}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{crit.fun}}, \code{\link{opt.par}}, \code{\link{opt.random.par}}, \code{\link{plot.mat}}}
\examples{
#generating a simple network corresponding to the simple Sum of squares
#structural equivalence with blockmodel:
# null com
# null null
n<-c(7,13)
net<-matrix(NA,nrow=n[1],ncol=n[2])
clu<-list(rep(1:2,times=c(3,4)),rep(1:2,times=c(5,8)))
tclu<-lapply(clu,table)
net[clu[[1]]==1,clu[[2]]==1]<-rnorm(n=tclu[[1]][1]*tclu[[2]][1],
   mean=0,sd=1)
net[clu[[1]]==1,clu[[2]]==2]<-rnorm(n=tclu[[1]][1]*tclu[[2]][2],
   mean=4,sd=1)
net[clu[[1]]==2,clu[[2]]==1]<-rnorm(n=tclu[[1]][2]*tclu[[2]][1],
   mean=4,sd=1)
net[clu[[1]]==2,clu[[2]]==2]<-rnorm(n=tclu[[1]][2]*tclu[[2]][2],
   mean=0,sd=1)
plot.mat(net,clu=clu) #two mode matrix of a two mode network
#converting to one mode network
M1<-two2one(net)$M
plot.mat(M1,clu=two2one(net)$clu) #plotting one mode matrix
plot.mat(one2two(M1,clu=clu)$M,clu=clu)
#converting one to two mode matix and ploting
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS

