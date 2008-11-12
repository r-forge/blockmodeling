\name{sedist}
\alias{sedist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes distances in terms of Structural equivalence (Lorrain and White, 1971)}
\description{
The functions computed the distances in terms of Structural equivalence (Lorrain and White, 1971) between the units of a one-mode network. Several options for treating the diagonal values are supported.
}
\usage{
sedist(M, method = "default", fun = "default", fun.on.rows = "default", handle.interaction = "switch", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network must be one-mode.}
  \item{method}{The method used to compute distances - any of the methods alloed by functions dist, cor or cov {all package::stats} or just "cor" or "cov" (given as character).}
  \item{fun}{Which function should be used to comput distacnes (given as character), .}
  \item{fun.on.rows}{For non-standard function - does the function compute measure on rows (such as \code{cor}, \code{cov},...) of the data matrix (as opposed to computing measure on columns (such as \code{dist}).}
  \item{handle.interaction}{How should the interaction between the vertices analysed be handled:\cr
						"switch" (the default) - assumes that when comparing units i and j, M[i,i] should be compared with M[j,j] and M[i,j] with M[j,i]\cr
						"ignore" (diagonal) - Diagonal is ignored\cr
						"none" - the matrix is used "as is"}
  \item{\dots}{Additional arguments to \code{fun}}
}
\details{
If both \code{method} and \code{fun} are "default", the euclidian distances are computed. the "default" method for \code{fun="dist"} is "eucludian" and for \code{fun="cor"} "pearson".
}
\value{
 A matrix (usually of class dist) is returned.
}
\references{
Batagelj, V., Ferligoj, A., Doreian, P. (1992): Direct and indirect methods for structural equivalence. Social Networks 14, 63-90.

Lorrain, F., White, H.C., 1971. Structural equivalence of individuals in social networks. Journal of Mathematical Sociology 1, 49-80.
}
\author{Aleš Žiberna}
\seealso{\code{\link{hist}}, \code{\link{REGE}}, \code{\link{crit.fun}},\code{\link{opt.par}},\code{\link{opt.random.par}},\code{\link{opt.these.par}},\code{\link{check.these.par}},\code{\link{plot.mat}}}
\examples{
#generating a simple network corresponding to the simple Sum of squares structural equivalence with blockmodel:
# null com
# null null
n<-20
net<-matrix(NA,ncol=n,nrow=n)
clu<-rep(1:2,times=c(5,15))
tclu<-table(clu)
net[clu==1,clu==1]<-rnorm(n=tclu[1]*tclu[1],mean=0,sd=1)
net[clu==1,clu==2]<-rnorm(n=tclu[1]*tclu[2],mean=4,sd=1)
net[clu==2,clu==1]<-rnorm(n=tclu[2]*tclu[1],mean=0,sd=1)
net[clu==2,clu==2]<-rnorm(n=tclu[2]*tclu[2],mean=0,sd=1)

D<-sedist(M=net)
plot.mat(net, clu=cutree(hclust(d=D,method="ward"),k=2))
}

\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS

