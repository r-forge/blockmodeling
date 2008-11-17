\name{REGE}
\alias{REGE}
\alias{REGE.for}
\alias{REGE.nm.for}
\alias{REGE.ow}
\alias{REGE.ow.for}
\alias{REGE.ownm.for}
\alias{REGD.for}
\alias{REGD.ow.for}
\alias{REGE.FC}
\alias{REGE.FC.ow}
\alias{REGD.ne.for}
\alias{REGD.ow.ne.for}
\alias{REGE.ne.for}
\alias{REGE.nm.diag.for}
\alias{REGE.nm.ne.for}
\alias{REGE.ow.ne.for}
\alias{REGE.ownm.diag.for}
\alias{REGE.ownm.ne.for}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{REGE - Algorithms for compiting (dis)similarities in terms of regular equivalnece.}
\description{
REGE - Algorithms for compiting (dis)similarities in terms of regular equivalnece (White and Reitz, 1983):

\code{REGE, REGE.for} - Classical REGE or REGGE, as also implemented in Ucinet. Similarities in terms of regular equivalnece are computed.  The \code{REGE.for} is a wraper for calling the FORTRAN subrutine writen by White (1985a), modified to be called by R. The \code{REGE} does the same, however it is written in R. The functions with and without ".for" differ only in wheater they are implemeted in R of FORTRAN. Needless to say, the funcitons implemeted in FORTRAN are much faster.

\code{REGE.ow, REGE.ow.for} - The above function, modified so that a best match is searhed for for each arc speleratly (and not for both arcs, if they exist, together)

\code{REGE.nm.for} - REGE or REGGE, modified to to use row and column normalited matrices instead of the original matrix.

\code{REGE.ownm.for} - The above function, modified so that a best match is searhed for for each arc speleratly (and not for both arcs, if they exist, together)

\code{REGD.for} - REGD or REGDI, a dissimilarity version of the classical REGE or REGGE. Dissimilarities in terms of regular equivalnece are computed.  The \code{REGD.for} is a wraper for calling the FORTRAN subrutine writen by White (1985b), modified to be called by R.

\code{REGE.FC}  - Acctually an erlier version of REGE. The diference is in the denominator. See Žiberna (2006) for details.

\code{REGE.FC.ow} - The above function, modified so that a best match is searhed for for each arc speleratly (and not for both arcs, if they exist, together)

other - still in testing stage
}
\usage{
REGE(M, E = 1, iter = 3, until.change = TRUE, use.diag = TRUE)
REGE.for(M, iter = 3, E = 1)
REGE.nm.for(M, iter = 3, E = 1)
REGE.ow(M, E = 1, iter = 3, until.change = TRUE, use.diag = TRUE)
REGE.ow.for(M, iter = 3, E = 1)
REGE.ownm.for(M, iter = 3, E = 1)
REGD.for(M, iter = 3, E = 0)
REGD.ow.for(M, iter = 3, E = 0)
REGE.FC(M, E = 1, iter = 3, until.change = TRUE, use.diag = TRUE,
   normE = FALSE)
REGE.FC.ow(M, E = 1, iter = 3, until.change = TRUE,
   use.diag = TRUE, normE = FALSE)
REGD.ne.for(M, iter = 3, E=0)
REGD.ow.ne.for(M, iter = 3, E = 0)
REGE.ne.for(M, iter = 3, E=1)
REGE.nm.diag.for(M, iter = 3, E=1)
REGE.nm.ne.for(M, iter = 3, E=1)
REGE.ow.ne.for(M, iter = 3, E=1)
REGE.ownm.diag.for(M, iter = 3, E=1)
REGE.ownm.ne.for(M, iter = 3, E=1)
}%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{Matrix or a 3 dimensional array representing the network. The third dimension allows for several relations to be analyzed.}
  \item{E}{Initial (dis)similarity in terms of regular equivalnece.}
  \item{iter}{The desired number of itetations}
  \item{until.change}{Should the iterations be stop when no change occours}
  \item{use.diag}{Should the diagonal be used. If \code{FALSE}, all diagonal elements are set to 0.}
  \item{normE}{Should the equivalence matrix be normalized after each iteration?}
}
\value{
  \item{E}{A matrix of (dis)similarities in terms of regular equivalnece}
  \item{Eall }{An array of (dis)similarity matrices in terms of regular equivalnece, each third dimmension represets one iteration. For ".for" functions, only the initial and the final (dis)similarities are returned.}
  \item{M}{Matrix or a 3 dimensional array representing the network used in the call.}
  \item{iter}{The desired number of itetations}
  \item{use.diag}{Should the diagonal be used - for functions implemeted in R only.}
  ...
}
\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

White, D. R., K. P. Reitz (1983): "Graph and semigroup homomorphisms on networks of relations". Social Networks, 5, p. 193-234.

White, Douglas R.(1985a): DOUG WHITE'S REGULAR EQUIVALENCE PROGRAM. \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/REGGE.FOR} (12.5.2005).

White, Douglas R.(1985b): DOUG WHITE'S REGULAR DISTANCES PROGRAM. \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/REGDI.FOR} (12.5.2005).

White, Douglas R.(2005): REGGE (web page). \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/} (12.5.2005).
}
\author{\enc{Aleš Žiberna}{Ales Ziberna} based on Douglas R. White's original REGE and REGD}
\seealso{\code{\link{sedist}}, \code{\link{crit.fun}}, \code{\link{opt.par}}, \code{\link{plot.mat}}}
\examples{
n<-20
net<-matrix(NA,ncol=n,nrow=n)
clu<-rep(1:2,times=c(5,15))
tclu<-table(clu)
net[clu==1,clu==1]<-0
net[clu==1,clu==2]<-rnorm(n=tclu[1]*tclu[2],
   mean=4,sd=1)*sample(c(0,1),
   size= tclu[1]*tclu[2],replace=TRUE,prob=c(3/5,2/5))
net[clu==2,clu==1]<-0
net[clu==2,clu==2]<-0

D<-REGE.for(M=net)$E #any other REGE function can be used
plot.mat(net, clu=cutree(hclust(d=as.dist(1-D),method="ward"),
   k=2))
#REGE returns similarities, which have to be converted to
#disimilarities

res<-opt.random.par(M=net,k=2,rep=10,approach="ss",blocks="reg",
   FUN="max")
plot(res) #Hopefully we get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
