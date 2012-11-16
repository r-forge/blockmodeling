\name{genRandomPar}
\alias{genRandomPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The function for generating random partitions}
\description{
 The function generates random partitions. The function is meant to be caled by the function \code{\link{opt.random.par}}
}
\usage{
genRandomPar(k, n, seed = NULL, mingr = 1, maxgr = Inf,
   addParam = list(genPajekPar = TRUE, probGenMech = NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{Number of clusters (by modes)}
  \item{n}{Number of units (by modes)}
  \item{seed}{Seed for generating random numbers (partitions)}
  \item{mingr}{Minimal alowed group size}
  \item{maxgr}{Maximal alowed group size}
  \item{addParam}{This has to be a list with the following parameters (any or all can be missing, then the default values (see usage) are used):\cr
  "genPajekPar" - Should the partitions be generated as in Pajek (Batagelj and Mrvar, 2006). If \code{FALSE}, all partitions are selected completly at random while making sure that the partitions have the required number of clusters. \cr
  "probGenMech" - Here the probabilities for 4 diferent generating mechanichems can be specified. If this is not specified, the value is set to \code{c(1/3,1/3,1/3,0)} if \code{genPajekPar} is \code{TRUE} and to \code{c(0,0,0,1)} if \code{genPajekPar} is \code{FALSE}. The first 3 mechanisems are the same as implemetned in Pajek (the second one has almost all units in only one cluster)  and the fourth is completly random (from uniform distribution).
}
}
\value{
  A random partition in the format required by \code{\link{opt.random.par}}. If a netowork has several modes, than a list of partitions, one for each mode.
}
\references{
BATAGELJ, Vladimir, MRVAR, Andrej (2006): Pajek 1.11, \url{http://vlado.fmf.uni lj.si/pub/networks/pajek/} (accessed January 6, 2006).
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{opt.random.par}}}
\keyword{cluster}
