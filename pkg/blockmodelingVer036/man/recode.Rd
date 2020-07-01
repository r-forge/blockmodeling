\name{recode}
\alias{recode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recode}
\description{
  Recodes values in a vector.
}
\usage{
recode(x, oldcode = sort(unique(x)), newcode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector.}
  \item{oldcode}{A vector of old codes.}
  \item{newcode}{A vector of new codes.}
}

\value{
  A recoded vector.
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\examples{
x <- rep(1:3, times = 1:3)
newx <- recode(x, oldcode = 1:3, newcode = c("a", "b", "c"))
}

\keyword{manip}% at least one, from doc/KEYWORDS

