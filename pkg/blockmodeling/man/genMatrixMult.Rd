\name{genMatrixMult}
\alias{genMatrixMult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized matrix multiplication}
\description{
Computes a generalized matrix multiplication, where sum and product functions (elemetwise and summary functions) can be replaced by arbitrary functions.}
\usage{
genMatrixMult(A,B,FUNelement="*", FUNsummary=sum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{The first matrix}
  \item{B}{The second matrix}
  \item{digits}{The number of desired digits}
  \item{FUN}{Function used for "shortening" the numbers.}
  \item{\dots}{Additional arguments to \code{format}}
}
\value{
A character vector or matrix.
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{matmult}}}
\examples{
#operations can be anything
x<-matrix(letters[1:8],ncol=2)
y<-matrix(1:10,nrow=2)

genMatrixMult(x,y,FUNelement=paste,FUNsummary=function(x)paste(x,collapse="|"))

#biary logic
set.seed(1)
x<-matrix(rbinom(8,size=1,prob=0.5)==1,ncol=2)
y<-matrix(rbinom(10,size=1,prob=0.5)==1,nrow=2)
genMatrixMult(x,y,FUNelement="*",FUNsummary=any)
}
\keyword{array} % at least one, from doc/KEYWORDS
\keyword{algebra}

