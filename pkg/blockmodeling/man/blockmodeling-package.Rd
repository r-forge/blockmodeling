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
Doreian, P., Batagelj, V. & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press.

\enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002

\enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207

\enc{Žiberna, A.}{Ziberna, A.} (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002

White, D. R., & Reitz, K. P. (1983). Graph and semigroup homomorphisms on networks of relations. Social Networks, 5(2), 193-234.

White, D. R. (2005). REGGE. Retrieved from http://eclectic.ss.uci.edu/~drwhite/REGGE/.
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
n <- 8 # If larger, the number of partitions increases dramatically,
# as does if we increase the number of clusters
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(3, 5))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

# We select a random partition  and then optimize  it
all.par <- nkpartitions(n = n, k = length(tclu))
# Forming the partitions
all.par <- lapply(apply(all.par, 1, list), function(x) x[[1]])

# Optimizing one partition
res <- optParC(M = net,
   clu = all.par[[sample(1:length(all.par), size = 1)]],
   approaches = "hom", homFun = "ss" , blocks = "com")
plot(res) # Hopefully we get the original partition

# Optimizing 10 random partitions with optRandomParC
res <- optRandomParC(M = net, k = 2, rep = 10, 
approaches = "hom", homFun = "ss", blocks = "com")
plot(res) # Hopefully we get the original partition

# Using indirect approach - structural equivalence
D <- sedist(M = net)
plot.mat(net, clu = cutree(hclust(d = D, method = "ward"), k = 2))
}
