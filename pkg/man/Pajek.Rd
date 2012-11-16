\name{Pajek}
\alias{Pajek}
\alias{loadnetwork}
\alias{loadnetwork2}
\alias{savevector}
\alias{savenetwork}
\alias{savematrix}
\alias{loadmatrix}
\alias{loadvector}
\alias{loadvector2}
\alias{loadpajek}
\title{Functions for loading and writing Pajek files}
\description{
  Functions for reading/loading and writing Pajek files:
 
  \code{loadnetwork} - Loads a Pajek ".net" filename as a matrix. For now, only simple one and two-mode networks are supported (eg. only single relations, no time information).
 
  \code{loadnetwork2} - The same as above, but adopted to be called withih \code{loadpajek}
 
  \code{savenetwork} - Saves a matrix in to a  Pajek ".net" filename.
 
  \code{loadmatrix} - Loads a Pajek ".mat" filename as a matrix.
 
  \code{savematrix} - Saves a matrix in to a  Pajek ".mat" filename.
 
  \code{loadvector} - Loads a Pajek ".clu" filename as a vector.
 
  \code{loadvector2} - The same as above, but adopted to be called withih \code{loadpajek} - as a consequence not suited for reading clusters
 
  \code{savevector} - Saves a vector in to a  Pajek ".clu" filename.

  \code{loadpajek} - Loads a Pajek project filename (".paj") as a list with the following components: Networks, Partitions, Vectors and Clusters. Clusters and hierarchies are dissmised.

}

\usage{
loadnetwork(filename,useSparseMatrix=NULL,minN=50)
loadnetwork2(filename,useSparseMatrix=NULL,minN=50)
savenetwork(n, filename, twomode = "default", symetric = NULL)
loadmatrix(filename)
savematrix(n, filename, twomode = 1)
loadvector(filename)
loadvector2(filename)
savevector(v, filename)
loadpajek(filename)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{The name of the filename to be loaded or saved to or an open file object.}
  \item{useSparseMatrix}{Should a sparse matrix be use instead of the ordinary one? Sparse matices can only be used if package Matrix is installed. The default \code{NULL} uses sparsematrices for networks with more that \code{minN} vertices}
  \item{minN}{The minimal number of units in the network to use sparse matrices.}
  \item{n}{A matrix representing the network.}
  \item{twomode}{1 for one-mode networks and 2 for two-mode networks. Default sets the argument to 1 for square matrices and to 2 for others.}
  \item{symetric}{If true, only the lower part of the matrix is used and the values are interpreted as "Edges", not "Arcs".}
  \item{v}{A vector}
}
\value{
  NULL, a matrix or a vector (see Description)
}
\references{
Pajek (
V. Batagelj, A. Mrvar: Pajek - Program for Large Network Analysis. Home page \url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}.

W. de Nooy, A. Mrvar, V. Batagelj: Exploratory Social Network Analysis with Pajek, CUP, January 2005
}
\author{Vladimir Batagelj & Andrej Mrvar (most functions), \enc{Aleš Žiberna}{Ales Ziberna} (\code{loadnetwork}, \code{loadpajek} and modification of others)}
\seealso{\code{\link{plot.mat}}, \code{\link{crit.fun}}, \code{\link{opt.par}}, \code{\link{opt.random.par}}, \code{\link{opt.these.par}}, \code{\link{check.these.par}}}
\keyword{graphs}% at least one, from doc/KEYWORDS
\keyword{file}% at least one, from doc/KEYWORDS
