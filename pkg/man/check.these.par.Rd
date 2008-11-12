\name{check.these.par}
\alias{check.these.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the value of a criterion function for a given network and a set of partitions}
\description{
The function computes the value of a criterion function for a given network and a set of partitions for Generalized blockmodeling. (Žiberna, 2006) based on other parameters (see below and \code{\link{crit.fun}}).
}
\usage{
check.these.par(M, partitions, approach, return.err = TRUE, save.initial.param = TRUE, force.fun = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{partitions}{A list of partitions. Each unique value represents one cluster. If the nework is  one-mode, than this should be a vector, else a list of vectors, one for each mode.}
  \item{approach}{One of the approaches described in Žiberna (2006). Possible values are:\cr
  "bin" - binary blockmodeling,\cr
  "val" - valued blockmodeling,\cr
  "imp" - implicit blockmodeling,\cr
  "ss" - sum of squares homogenity blockmodeling, and\cr
  "ad" - absolute deviations homogenity blockmodeling.}
  \item{return.err}{Should the error for each evaluated partition be returned}
  \item{save.initial.param}{Should the inital parameters (\code{approach},...)}
  \item{force.fun}{Select the function used to evaluate the network and a partition. This should be used only in exterem cases. Otherwise, the appropriate function is selected (generated) besed on the input parameters.}
  \item{\dots}{Argumets to \code{gen.crit.fun} see \code{\link{crit.fun}} for description. Some are required!!!}
}
\value{
  \item{M}{The matrix of the network analyzed}
  \item{best}{A list of results from \code{crit.fun.tmp} with the same elements as the result of \code{crit.fun}, only without \code{M}}
  \item{err}{If selected - The vector of errors or inconsistencies of the emplirical network with the ideal network for a given blockmodel (model,approach,...) and parititions}
  \item{call}{The call used to call the function.}
  \item{initial.param}{If selected - The inital parameters used.}
    ...
}
\section{Warning }{
This function is usually used to check all possible partitions. If the number of partitions is large (several 1000), this can be extremly time demanding. It is advaisable to firtst time the function on a smaller subset.
}

\references{
ŽIBERNA, Aleš (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

ŽIBERNA, Aleš. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, Anuška (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.
}

\author{Aleš Žiberna}
\seealso{\code{\link{crit.fun}}, \code{\link{opt.par}}, \code{\link{opt.random.par}}, \code{\link{opt.these.par}}, \code{\link{nkpartitions}}, \code{\link{nkpar}}, \code{\link{plot.check.these.par}} }

\examples{
n<-8 #if larger, the number of partitions increases dramaticaly, as does if we increase the number of clusters
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
all.par<-lapply(apply(all.par,1,list),function(x)x[[1]]) # to make a list out of the matrix
res<-check.these.par(M=net,partitions=all.par,approach="ss",blocks="com")
plot(res) #we get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
