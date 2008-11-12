\name{rand}
\alias{crand}
\alias{crand2}
\alias{rand}
\alias{rand2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparing partitions}
\description{
Rand Index and Rand Index corrected/adjusted for chance for comparing partitions (Hubert and Arabie, 1985). The names of the clusters do not matter.
}
\usage{
rand(tab)
rand2(clu1, clu2)
crand(tab)
crand2(clu1, clu2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clu1, clu2}{The two partitions to be compared, given in the from of vectors, where for each unit a cluster membership is given.}
  \item{tab}{A contingency table obtaind as table(clu1,clu2)}
}
\value{
  The value of Rand Index (corrected/adjusted for chance)
}
\references{Hubert L. in Arabie P. (1985): Comparing Partitions. Journal of Classification, 2, 193-218.}
\author{Aleš Žiberna}
\keyword{cluster}% at least one, from doc/KEYWORDS
