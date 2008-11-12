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
\title{Functions for loading and writing Pajek[http://vlado.fmf.uni-lj.si/pub/networks/pajek/] files}
\description{
  Functions for reading/loading and writing Pajek files:
  
  \code{loadnetwork} - Loads a Pajek ".net" file as a matrix. For now, only simple one and two-mode networks are supported (eg. only single relations, no time information).
  
  \code{loadnetwork2} - The same as above, but adopted to be called withih \code{loadpajek}
  
  \code{savenetwork} - Saves a matrix in to a  Pajek ".net" file.
  
  \code{loadmatrix} - Loads a Pajek ".mat" file as a matrix.
  
  \code{savematrix} - Saves a matrix in to a  Pajek ".mat" file.
  
  \code{loadvector} - Loads a Pajek ".clu" file as a vector.
  
  \code{loadvector2} - The same as above, but adopted to be called withih \code{loadpajek} - as a consequence not suited for reading clusters
  
  \code{savevector} - Saves a vector in to a  Pajek ".clu" file.

  \code{loadpajek} - Loads a Pajek project file (".paj") as a list with the following components: Networks, Partitions, Vectors and Clusters. Clusters and hierarchies are dissmised.

}

\usage{
loadnetwork(file,useSparseMatrix=NULL,minN=50)
loadnetwork2(file,useSparseMatrix=NULL,minN=50)
savenetwork(n, file, twomode = "default", symetric = all(n==t(n)))
loadmatrix(file)
savematrix(n, file, twomode = 1)
loadvector(file)
loadvector2(file)
savevector(v, file)
loadpajek(file)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{The name of the file to be loaded or saved to or an open file object.}
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
V. Batagelj, A. Mrvar: Pajek - Program for Large Network Analysis. Home page
http://vlado.fmf.uni-lj.si/pub/networks/pajek/.

W. de Nooy, A. Mrvar, V. Batagelj: Exploratory Social Network Analysis with Pajek, CUP, January 2005
}
\author{Vladimir Batagelj & Andrej Mrvar (most functions), Aleš Žiberna (\code{loadnetwork},\code{loadpajek} and modification of others)}
\seealso{\code{\link{plot.mat}},\code{\link{crit.fun}},\code{\link{opt.par}},\code{\link{opt.random.par}},\code{\link{opt.these.par}},\code{\link{check.these.par}}}
\keyword{graphs}% at least one, from doc/KEYWORDS
\keyword{file}% at least one, from doc/KEYWORDS
