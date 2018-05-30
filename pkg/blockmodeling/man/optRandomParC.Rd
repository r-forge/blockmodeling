
\name{optRandomParC}
\alias{optRandomParC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimizing a set of partitions based on the value of a criterion function}
\description{
The function optimizes a set of partitions based on the value of a criterion function  (see \code{\link{critFunC}} for details on the criterion function) for a given network and blockmodel for Generalized blockmodeling (Žiberna, 2007) based on other parameters (see below).
The optimization is done through local optimization, where the neighborhood of a partition includes all partitions that can be obtained by moving one unit from one cluster to another or by exchanging two units (from different clusters).
A list of paritions can or the number of clusters and a number of partitions to generate can be specified (\code{optParC}).

}
\usage{

optRandomParC(M, k, rep, save.initial.param = TRUE,
              save.initial.param.opt = FALSE, deleteMs = TRUE,
              max.iden = 10, switch.names = NULL, return.all =
              FALSE, return.err = TRUE, seed = NULL, RandomSeed =
              NULL, parGenFun = genRandomPar, mingr = NULL, maxgr =
              NULL, addParam = list(genPajekPar = TRUE, probGenMech
              = NULL), maxTriesToFindNewPar = rep * 10, skip.par =
              NULL, useOptParMultiC = FALSE, useMulti =
              useOptParMultiC, printRep = ifelse(rep <= 10, 1,
              round(rep/10)), n = NULL, nCores = 1, ...) 


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks   are supported. The network can have one or more modes (different kinds of units with no ties among themselves. If the network is not two-mode, the matrix must be square.}
  \item{k}{The number of clusters used in the generation of partitions.}
  \item{rep}{The number of repetitions/different starting partitions to check.}
  \item{save.initial.param}{Should the initial parameters (\code{approaches}, ...) be saved. The        default value is \code{TRUE}.}
  \item{save.initial.param.opt}{Should the inital parameters(\code{approaches}, ...) of using          \code{optParC} be saved. The default value is \code{FALSE}.}
  \item{deleteMs}{Delete networks/matrices from the results of to save space.}
  \item{max.iden}{Maximum number of results that should be saved (in case there are more than          \code{max.iden} results with minimal error, only the first \code{max.iden} will be saved).}
  \item{switch.names}{Should partitions that only differ in group names be considered equal.}
  \item{return.all}{If \code{FALSE}, solution for only the best (one or more) partition/s is/are returned.}
  \item{return.err}{Should the error for each optimized partition be returned.}
  \item{seed}{Optional. The seed for random generation of partitions.}
  \item{RandomSeed}{Optional. Integer vector, containing the random number generator.It is only         looked for in the user's workspace.}
  \item{parGenFun}{The function (object) that will generate random partitions. The default function is   \code{\link{genRandomPar}}. The function has to accept the following parameters: \code{k} (number o of partitions by modes, \code{n} (number of units by modes), \code{seed} (seed value for random generation of partition), \code{addParam} (a list of additional parameters).}
  \item{mingr}{Minimal allowed group size.}
  \item{maxgr}{Maximal allowed group size.}
  \item{addParam}{A list of additional parameters for function specified above. In the usage section they are specified for the default function \code{\link{genRandomPar}}.}
  \item{genPajekPar}{Should the partitions be generated as in Pajek.}
  \item{probGenMech}{Should the probabilities for different mechanisms for specifying the partitions be set. If \code{probGenMech} is not set, it is determined based on the parameter \code{genPajekPar}.}
  \item{maxTriesToFindNewPar}{The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before - the default value is \code{rep*1000}.}
  \item{skip.par}{The partitions that are not allowed or were already checked and should therefore be skipped.}
  \item{useOptParMultiC}{For backward compatibility. May be removed soon. See next argument.}
  \item{useMulti}{Which version of local search should be used. Default is currently \code{FALSE}. If \code{FALSE}, first possible all moves in random order and then all possible exchanges in random order are tired. When a move with lower value of criterion function is found, the algorithm moves to this new partition. If \code{TRUE} the version of local search where all possible moves and exchanges are tired first and then the one with the lowest error is selected and used. In this case, several optimal partitions are found. \code{maxPar} best partitions are returned.}
 \item{printRep}{Should some information about each optimization be printed.}
 \item{n}{The number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is one mode (both modes are in rows and columns).}
  \item{nCores}{Number of cores to be used. Value \code{0} means all available cores. It can also be a cluster object.}
  \item{\dots}{Arguments passed to other functions, see \code{\link{critFunC}}.}
}

\value{
  \item{M}{The matrix of the network analyzed.}
  \item{res}{If \code{return.all = TRUE} - A list of results the same as \code{best} - one \code{best} for each partition optimized.}
  \item{best}{A list of results from \code{crit.fun.tmp} with the same elements as the result of \code{crit.fun}, only without \code{M}.}
  \item{err}{If \code{return.err = TRUE} - The vector of errors or inconsistencies of the empirical  network with the ideal network for a given blockmodel (model,approach,...) and parititions.}
  \item{nIter}{The vector of the number of iterations used - one value for each starting partition that was optimized. It can show that \code{maxiter} is too low if a lot of these values have the value of \code{maxiter}.}
  \item{checked.par}{If selected - A list of checked partitions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includes the partitions in \code{skip.par}.}
  \item{call}{The call used to call the function.}
  \item{initial.param}{If selected - The initial parameters are used.}
}
\section{Warning}{
It should be noted that the time complexity of package blockmodeling is increasing with
the number of units and the number of clusters (due to its algorithm). Therefore the analysis
of network with more than 100 units can take a lot of time (from a few hours to a few days).
}

\references{
Batagelj, V., & Mrvar, A. (2006). Pajek 1.11. Retrieved from \url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}

Doreian, P., Batagelj, V. & Ferligoj, A. \enc{Anuška}{Anuska}. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press.

\enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002

\enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207

\enc{Žiberna, A.}{Ziberna, A.} (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002

}

\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{critFunC}}}
\examples{
n <- 8 # If larger, the number of partitions increases dramatically 
# as does if we increase the number of clusters
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(3, 5))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

# We select a random parition and then optimise it
all.par <- nkpartitions(n = n, k = length(tclu))
# Forming the partitions
all.par <- lapply(apply(all.par, 1, list), function(x)x[[1]])

# Optimizing one partition
res <- optParC(M = net,
   clu = all.par[[sample(1:length(all.par), size = 1)]],
   approaches = "hom", homFun = "ss", blocks = "com")
plot(res) # Hopefully we get the original partition

# Optimizing 10 random chosen partitions with optRandomParC
res <- optRandomParC(M = net, k = 2, rep = 10, 
approaches = "hom", homFun = "ss", blocks = "com")
plot(res) # Hopefully we get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
