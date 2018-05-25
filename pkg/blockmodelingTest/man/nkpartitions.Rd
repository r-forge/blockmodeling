\name{nkpartitions}
\alias{nkpartitions}
\alias{nkpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for listing all possible partitions or just counting  the number of them}
\description{
The function \code{nkpartitions} lists all possible partitions of n objects in to k clusters. The function \code{nkpar} only gives the number of such partitions.
}
\usage{
nkpartitions(n, k, exact = TRUE, print = FALSE)
nkpar(n, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of units/objects.}
  \item{k}{Number of clusters/groups.}
  \item{exact}{Search for partitions with exactly \code{k} or at most \code{k} clusters.}
  \item{print}{Print results as they are found.}
}
\value{
The matrix or number of possible partitions.
}
\author{Chris Andrews}
\examples{
n <- 8 # If larger, the number of partitions increases dramaticaly,
     # as does if we increase the number of clusters
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(3, 5))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

# Computation of criterion function with the correct partition
nkpar(n = n, k = length(tclu)) # Computing the number of partitions
all.par <- nkpartitions(n = n, k = length(tclu)) # Forming the partitions
all.par <- lapply(apply(all.par, 1, list), function(x) x[[1]])
# to make a list out of the matrix
res <- critFunC(M = net, clu = clu, approaches = "val", 
blocks = c("nul", "com"), preSpecM = 4)
plot(res) # We get the original partition
}
\keyword{cluster}% at least one, from doc/KEYWORDS

