\name{ircNorm}
\alias{ircNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for iterated row and column normalization of valued matrices.}
\description{
 The aim is to obtain a matrix with row and column sums equal to 1. This is achieved by iterating row and column normalization. This is usually not possible if any row or column has only 1 non-zero cell.
}
\usage{
ircNorm(M, eps = 10^-12, maxiter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A non-negative valued matrix to be normalized}
  \item{eps}{The maximum allows squared deviation of a row or column maximum from 1 (if not exaclty 0). Also, if the all deviations in  to consequtive iterations are smaller, the process is termianted.}
  \item{maxiter}{Maximum number of iterations. If reached the process is termianted and the current solution returned}
}
\value{
  Normalized matrix.
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\examples{

A<-matrix(runif(100),ncol=10)
A #A non-normalized matrix with different row and column sums.
apply(A,1,sum)
apply(A,2,sum)
A.norm<-ircNorm(A)
A.norm #Normalized matrix with all row and column sums aproximately 1.
apply(A.norm,1,sum)
apply(A.norm,2,sum)
}

\keyword{manip}

