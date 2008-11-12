\name{opt.random.par}
\alias{opt.random.par}
\alias{opt.these.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimizes a set of partitions based on the value of a criterion function.}
\description{
The function optimizes a set partitions based on the value of a criterion function  (see \code{\link{crit.fun}} for details on the criterion function) for a given network and blockmodel for Generalized blockmodeling (Žiberna, 2006) based on other parameters (see below).
The optimization is done through local optimization, where the neighborhood of a partition includes all partitions that can be obtained by moving one unit from one cluster to another or by exchanging two units (from different clusters).
A list of paritions can be specified (\code{opt.these.par}) or the number of clusters and a number of partitions to generate (\code{opt.random.par}).
}\usage{
opt.random.par(M, k, n = NULL, rep, return.all = FALSE, return.err = TRUE, maxiter = 50, m = NULL, approach, trace.iter = FALSE, switch.names = NULL, save.initial.param = TRUE, skip.par = NULL, save.checked.par = TRUE, merge.save.skip.par = any(!is.null(skip.par), save.checked.par), skip.allready.checked.par = TRUE, check.skip = "iter", print.iter = FALSE, max.iden = 10, seed = NULL, parGenFun = genRandomPar, mingr = 1, maxgr = Inf, addParam = list(genPajekPar = TRUE, probGenMech = NULL), maxTriesToFindNewPar = rep * 10, ...)

opt.these.par(M, partitions, return.all = FALSE, return.err = TRUE, skip.allready.checked.par = TRUE, maxiter = 50, m = NULL, approach, trace.iter = FALSE, switch.names = TRUE, save.initial.param = TRUE, skip.par = NULL, save.checked.par = !is.null(skip.par), merge.save.skip.par = all(!is.null(skip.par), save.checked.par), check.skip = "never", print.iter = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{k}{The number of clustrs used in generation of partitions.}
  \item{n}{The vector of the number of units in each mode (only necessary if mode is larger than 2.}
  \item{rep}{The number of repetitions/different starting partitions to check.}
  \item{partitions}{A list of partitions. Each unique value represents one cluster. If the nework is  one-mode, than this should be a vector, else a list of vectors, one for each mode.}
  \item{return.all}{If \code{FALSE}, solution for only the best (one or more) partition/s is/are returned.}
  \item{return.err}{Should the error for each optimized partition be returned}
  \item{maxiter}{Maximum number of iterations}
  \item{m}{Suficient value for individual cells for valued approach. Can be a number or a character string giving the name of a function. Set to \code{"max"} for implicit approach.}
  \item{approach}{One of the approaches described in Žiberna (2007). Possible values are:\cr
  "bin" - binary blockmodeling,\cr
  "val" - valued blockmodeling,\cr
  "imp" - implicit blockmodeling,\cr
  "ss" - sum of squares homogenity blockmodeling, and\cr
  "ad" - absolute deviations homogenity blockmodeling.}
  \item{trace.iter}{Should the result of each iteration (and not only of the best one) be saved}
  \item{switch.names}{Should  partitions that differ only in diferent names of positions be treated as different. It should be set to \code{TRUE} only if a asymetric blockmodel via \code{BLOCKS} is specified.}
  \item{save.initial.param}{Should the inital parameters (\code{approach},...) be saved}
  \item{skip.par}{The partitions that are not allowed or were already checked and should therfire be skiped.}
  \item{save.checked.par}{Should the checked partitions be saved. For example, so that they can be used in the next call as \code{skip.par}}
  \item{merge.save.skip.par}{Should the checked partitions be merged with skiped ones?}
  \item{skip.allready.checked.par}{If \code{TRUE},the partitions that were already checked when runing \code{opt.par} form different statrting points will be skiped.}
  \item{check.skip}{When should the check be preformed:\cr
  "all"  - before every call to 'crit.fun' (Time demanding)\cr
	"iter" - at the end of eack iteratiton\cr
  "opt.par"  - before every call to \code{opt.par}, when starting the optimization of a new partition.\cr
  "never" - never}
  \item{print.iter}{Should the progress of each iteration be printed?}
  \item{max.iden}{The maximum number of results that should be saved (in case there are more than max.iden results with minimal error, only the first max.iden will be saved).}
  \item{seed}{Optional. The seed for random generation of partitions.}
  \item{parGenFun}{The fucntion (object) that will generate rendom partitions. The deault fuction is \code{\link{genRandomPar}}. The function has to accept the following parameters: k (number of partitions by modes, n (number of units by modes), seed (seed value for random generation of partition), addParam (a list of additional parametres)}
  \item{mingr}{Minimal alowed group size}
  \item{maxgr}{Maximal alowed group size}
  \item{addParam}{A list of additional parameters for function specified above. In the usage section they are specified for the dthe default function \code{\link{genRandomPar}}:}
    \item{maxTriesToFindNewPar}{The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before - the default value is \code{rep*1000}}
  \item{\dots}{Argumets passed to other functions, see \code{\link{crit.fun}}}
}

\value{
  \item{M}{The matrix of the network analyzed}
  \item{res}{If \code{return.all = TRUE} - A list of results the same as \code{best} - one \code{best} for each partition optimized}
  \item{best}{A list of results from \code{crit.fun.tmp} with the same elements as the result of \code{crit.fun}, only without \code{M}}
  \item{err}{If \code{return.err = TRUE} - The vector of errors or inconsistencies of the emplirical network with the ideal network for a given blockmodel (model,approach,...) and parititions}
  \item{nIter}{The vector of number of iterations used - one value for each starting partition that was optimized. It can show that \code{maxiter} is to low if a lot of these values have the value of \code{maxiter}}
  \item{checked.par}{If selected - A list of checked parititions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includs the partitions in \code{skip.par}.}
  \item{call}{The call used to call the function.}
  \item{initial.param}{If selected - The inital parameters used.}
}
\section{Warning }{
This function can be extremly slow. The time complexity is incrising with the number od units and the number of clusters. It is advaisable to firtst time the function on a smaller network.
}

\references{
ŽIBERNA, Aleš (2007): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

ŽIBERNA, Aleš. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, Anuška (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.

BATAGELJ, Vladimir, MRVAR, Andrej (2006): Pajek 1.11, http://vlado.fmf.uni lj.si/pub/networks/pajek/ (accessed January 6, 2006).
}

\author{Aleš Žiberna}
\seealso{\code{\link{crit.fun}},\code{\link{check.these.par}},\code{\link{opt.par}},\code{\link{plot.opt.more.par}}}
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

#optimizing one partition
res<-opt.par(M=net,clu=all.par[sample(1:length(all.par),size=1)],approach="ss",blocks="com")
plot(res) #Hopefully we get the original partition

#optimizing 10 random chosen partitions which are suplied to opt.these.par
res<-opt.these.par(M=net,partitions=all.par[sample(1:length(all.par),size=10)],approach="ss",blocks="com")
plot(res) #Hopefully we get the original partition

#optimizing 10 random chosen partitions with opt.random.par
res<-opt.random.par(M=net,k=2,npar=10,approach="ss",blocks="com")
plot(res) #Hopefully we get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
