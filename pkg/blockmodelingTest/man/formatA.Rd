\name{formatA}
\alias{formatA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A formating function for numbers}
\description{
Formats a vector or matrix of numbers so that all have equal length (digits). This is especially suitable for printing tables.
}
\usage{
formatA(x, digits = 2, FUN = round, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numerical vector or matrix.}
  \item{digits}{The number of desired digits.}
  \item{FUN}{Function used for "shortening" the numbers.}
  \item{\dots}{Additional arguments to \code{format}.}
}
\value{
A character vector or matrix.
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{find.m}}, \code{\link{find.m2}}, \code{\link{find.cut}}}
\examples{
A <- matrix(c(1, 1.02002, 0.2, 10.3), ncol = 2)
formatA(A)
}
\keyword{character}% at least one, from doc/KEYWORDS
