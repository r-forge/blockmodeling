\name{clu}
\alias{clu}
\alias{partitions}
\alias{IM}
\alias{err}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for extraction of some elements for objects, returend by functions for Generalized blockmodeling}
\description{
Function for extraction of clu (partition), all best clus (partitions), IM (image or blockmodel) and err (total error or inconsistency) for objects, returend by functions \code{\link{critFunC}}or \code{\link{optRandomParC}}}

\usage{
clu(res, which = 1, ...)
IM(res, which = 1, ...)
err(res, ...)
partitions(res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{Result of function \code{\link{critFunC}} or \code{\link{optRandomParC}}}
  \item{which}{From which (if there are more than one) "best" solution whould the element be extracted. Warning! \code{which} grater than the number of "best" partitions produces an error.}
  \item{\dots}{Not used}
}
\value{
  The desired element.
}

\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, \enc{Anuška}{Anuska} (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.
}

\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{critFunC}}, \code{\link{plot.mat}}, \code{\link{optRandomParC}}}
\examples{
n <- 8 # if larger, the number of partitions increases dramaticaly,
     # as does if we increase the number of clusters
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(3, 5))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

# we select a random parition and then optimise it

all.par <- nkpartitions(n = n, k = length(tclu))
# forming the partitions
all.par <- lapply(apply(all.par, 1, list),function(x) x[[1]])
# to make a list out of the matrix
res <- optParC(M = net,
   clu = all.par[[sample(1:length(all.par), size = 1)]],
    approaches = "hom", homFun = "ss", blocks = "com")
plot(res) # Hopefully we get the original partition
clu(res) # Hopefully we get the original partition
err(res) # Error
IM(res) # Image matrix/array.
}
\keyword{manip}% at least one, from doc/KEYWORDS
