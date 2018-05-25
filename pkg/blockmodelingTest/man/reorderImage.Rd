\name{reorderImage}
\alias{reorderImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reordering an image matrix of the blockmodel (or an error matrix based on new and old partition}
\description{Reorders an image matrix of the blockmodel (or an error matrix based on new and old partition. The partitions should be the same, except that classes can have different labels. It is useful when we want to have a different order of classes in figures and then also in image matrices. Currently it is only suitable for one-mode blockmodels.
}
\usage{
reorderImage(IM, oldClu, newClu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IM}{An image or error matrix.}
  \item{oldClu}{Old partition.}
  \item{newClu}{New partition, the same as the old one except for class labeles.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Reorder matrix (rows and columns are reordred)
}
\references{
\enc{Žiberna, A.}.{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002

\enc{Žiberna, A.}.{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
}
\author{Ales Ziberna}
\seealso{\code{\link{critFunC}}, \code{\link{plot.mat}}, \code{\link{clu}}, \code{\link{IM}}, \code{\link{err}}}
%\examples{
%}
\keyword{manip}% at least one, from doc/KEYWORDS
