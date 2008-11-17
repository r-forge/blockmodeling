\name{opt.par}
\alias{opt.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimizes a partition based on the value of a criterion function.}
\description{
The function optimizes a partition based on the value of a criterion function (see \code{\link{crit.fun}}) for a given network and blockmodel for Generalized blockmodeling (Žiberna, 2006) based on other parameters (see below).
The optimization is done through local optimization, where the neighbourhood of a partition includes all partitions that can be obtained by moving one unit from one cluster to another or by exchanging two units (from different clusters).
}
\usage{
opt.par(M, clu, maxiter = 50, m = NULL, approach,
   trace.iter = FALSE, switch.names = NULL,
   save.initial.param = TRUE, skip.par = NULL,
   save.checked.par = !is.null(skip.par),
   merge.save.skip.par = all(!is.null(skip.par),
   save.checked.par), check.skip = "never", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{clu}{A partition. Each unique value represents one cluster. If the nework is  one-mode, than this should be a vector, else a list of vectors, one for each mode}
  \item{maxiter}{Maxsimum number of iterations allowed}
  \item{m}{Suficient value for individual cells for valued approach. Can be a number or a character string giving the name of a function. Set to \code{"max"} for implicit approach.}
  \item{approach}{One of the approaches described in Žiberna (2006). Possible values are:\cr
  "bin" - binary blockmodeling,\cr
  "val" - valued blockmodeling,\cr
  "imp" - implicit blockmodeling,\cr
  "ss" - sum of squares homogenity blockmodeling, and\cr
  "ad" - absolute deviations homogenity blockmodeling.}
  \item{trace.iter}{Should the result of each iteration (and not only of the best one) be saved}
  \item{switch.names}{Should  partitions that differ only in diferent names of positions be treated as different. It should be set to \code{TRUE} only if a asymetric blockmodel via \code{BLOCKS} is specified. The default \code{NULL} tries to find that.}
  \item{save.initial.param}{Should the inital parameters (\code{approach},...) be saved}
  \item{skip.par}{The partitions that are not allowed or were already checked and should therfire be skiped.}
  \item{save.checked.par}{Should the checked partitions be saved. For example, so that they can be used in the next call as \code{skip.par}}
  \item{merge.save.skip.par}{Should the checked partitions be merged with skiped ones?}
  \item{check.skip}{When should the check be preformed:\cr
  "all"  - before every call to 'crit.fun'\cr
  "iter" - at the end of eack iteratiton\cr
  "never" - never}
  \item{use.for.opt}{Should FORTRAN function be used for optimization if possible. If FORTRAN function is used, the speed is dramatically increast, however some the output is slightly different and the plotting function might not work. FORTRAN subrutines are available for only very special cases, currently only for using "ss" aproach and only complete blocks. If you are using such setting and some special features (these are not implemented in FORTRAN subrutines - e.g. using function \code{parOK} to allow only certain kinds of partitions), it's safer to set it to FASLE, as the fuction may miss that these features are not implemented in FORTRAN subrutines and use them nevertheless, leading to wrong results.}
    \item{use.for}{(default = \code{TRUE}) Should FORTRAN subrutines be used where available (available for only very special cases, currently only for using "ss" aproach and only complete blocks. If you are using such setting and some special features (these are not implemented in FORTRAN subrutines), it's safer to set it to FASLE, as the fuction may miss that these features are not implemented in FORTRAN subrutines and use them nevertheless, leading to wrong results.}
  \item{check.switch}{If \code{TRUE} (the default), the neighborhood of the selected partition also includes the partitions that can be obtained by exchanging (switching) two units from diferent clusters).}
  \item{check.all}{If \code{TRUE} (the default), all partitions in the neighborhood of the selected partition are first evaluated and the current partition than changes to the one with the lowest value of the criterion function (if lower than that of the current partition). If \code{FALSE}, the first partition with the criterion lower the current partition becomes the new current partition (and the iteration terminates).}
  \item{\dots}{Argumets passed to other functions, see \code{\link{crit.fun}}}
}
\value{
  \item{M}{The matrix of the network analyzed}
  \item{best}{A list of results from \code{crit.fun.tmp} with the same elements as the result of \code{crit.fun}, only without \code{M}}
  \item{iter}{A list of resoults the same as \code{best} - one \code{best} for each iteration}
  \item{err}{If selected - The vector of errors or inconsistencies of the emplirical network with the ideal network for a given blockmodel (model,approach,...) and parititions}
  \item{nIter}{The number of iterations used. It can show that \code{maxiter} is to low if this value is equal to \code{maxiter}}
  \item{call}{The call used to call the function.}
  \item{initial.param}{If selected - The inital parameters used.}
  \item{checked.par}{If selected - A list of checked parititions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includs the partitions in \code{skip.par}.}
    ...
}
\section{Warning }{
This function can be extremly slow. The time complexity is incrising with the number od units and the number of clusters. It is advaisable to firtst time the function on a smaller network.
}

\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, \enc{Anuška}{Anuska} (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.
}

\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{crit.fun}}, \code{\link{check.these.par}}, \code{\link{opt.random.par}}, \code{\link{opt.these.par}}, \code{\link{plot.opt.par}}}
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

#we select a random parition and then optimise it

all.par<-nkpartitions(n=n, k=length(tclu)) #forming the partitions
all.par<-lapply(apply(all.par,1,list),function(x)x[[1]])
# to make a list out of the matrix
res<-opt.par(M=net,
   clu=all.par[[sample(1:length(all.par),size=1)]],
   approach="ss",blocks="com")
plot(res) #Hopefully we get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
