\name{fun.by.blocks}
\alias{fun.by.blocks}
\alias{fun.by.blocks.default}
\alias{fun.by.blocks.mat}
\alias{fun.by.blocks.opt.more.par}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of function values by blocks}
\description{
Computes a value of a functions over blocks of a matrix, defined by a partition.
}
\usage{
fun.by.blocks(x, ...)

\method{fun.by.blocks}{default}(x = M, M = x, clu, 
ignore.diag = "default", sortNames = TRUE, 
FUN = "mean", ...) 

\method{fun.by.blocks}{opt.more.par}(x, which = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of suitable class or a matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{clu}{A partition. Each unique value represents one cluster. If the nework is  one-mode, than this should be a vector, else a list of vectors, one for each mode}
  \item{ignore.diag}{Should the diagonal be ingored. }
  \item{sortNames}{Should the rows and columns of the matrix be sorted based on their names?}
  \item{FUN}{Function to be computed over the blocks}
  \item{which}{Which (if several) of the "best" solutions should be used}
  \item{\dots}{Further arguments to \code{fun.by.blocks.default} }

}

\value{
 A numerical matrix of \code{FUN} values by blocks, induced by a partition \code{clu}
}
\references{
\enc{Žiberna, A.}.{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002

\enc{Žiberna, A.}.{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
}

\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{optRandomParC}}, \code{\link{optParC}}}

\examples{
n <- 8 #if larger, the number of partitions increases dramaticaly,
     #as does if we increase the number of clusters
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(3, 5))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)


#optimizing 10 random partitions with optRandomParC
res <- optRandomParC(M = net, k = 2, rep = 10, approaches = "hom", homFun = "ss", blocks = "com")
plot(res) #Hopefully we get the original partition
fun.by.blocks(res)
#computing mean by blocks, ignoring the diagonal (default)
}
\keyword{ cluster }% at least one, from doc/KEYWORDS
\keyword{ math }% at least one, from doc/KEYWORDS
