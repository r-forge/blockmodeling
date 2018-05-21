\name{blockmodeling-package}
\alias{blockmodeling-package}
\docType{package}
\title{
An R package for Generalized and classical blockmodeling of valued networks
}
\description{
This package is primarily meant as an implementation of Generalized blockmodeling. In addition, functions for computation of (dis)similarities in terms of structural and regular equivalence, plotting and other "utility" functions are provided.
}
%\details{
%~~ An overview of how to use the package, including the most important functions ~~
%}
\author{
\enc{Aleš Žiberna}{Ales Ziberna}
}
\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2014): Blockmodeling of multilevel networks. Social Networks, Oct. 2014, vol. 39, no. 1, 46-61. \url{https://doi.org/10.1016/j.socnet.2014.04.002}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, \enc{Anuška}{Anuska} (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.

White, D. R., K. P. Reitz (1983): "Graph and semigroup homomorphisms on networks of relations". Social Networks, 5, p. 193-234.

White, Douglas R.(2005): REGGE (web page). \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/} (12.5.2005).
}
\keyword{ package }
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS

\seealso{
Packages:
\code{\link[sna:sna]{sna}}
\code{\link[network:network-package]{network}}

Functions inside this package:
\code{\link{critFunC}}, \code{\link{optParC}}, \code{\link{optRandomParC}},  \code{\link{REGE}}, \code{\link{plot.mat}}
}
\examples{
n <- 8 #if larger, the number of partitions increases dramaticaly,
# as does if we increase the number of clusters
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(3, 5))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

#we select a random parition and then optimise it

all.par<-nkpartitions(n=n, k=length(tclu))
#forming the partitions
all.par <- lapply(apply(all.par, 1, list), function(x) x[[1]])
# to make a list out of the matrix

#optimizing one partition
res <- optParC(M = net,
   clu = all.par[[sample(1:length(all.par), size = 1)]],
   approaches = "hom", homFun = "ss" , blocks = "com")
plot(res) #Hopefully we get the original partition


#optimizing 10 random partitions with opt.random.par
res <- optRandomParC(M = net, k = 2, rep = 10, 
approaches = "hom", homFun = "ss", blocks = "com")
plot(res) #Hopefully we get the original partition



#using indidect approach - structural equivalence
D <- sedist(M = net)
plot.mat(net, clu = cutree(hclust(d = D, method = "ward"), k = 2))
}
