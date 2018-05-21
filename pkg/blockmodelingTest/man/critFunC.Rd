\name{critFunC}
\alias{critFunC}
\alias{optParC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for Generalized blockmodeling for valued networks.}
\description{
Functions for implementation of Generalized blockmodeling for valued
networks where the values of the ties are assumed to be measured on at least interval
scale. \code{critFunC} calculates criterion function, based on the network, partition and blockmodel/equivalece. \code{optParC} optimizes a partition based on the criterion function based on a local search algorithm.
}
\usage{
critFunC(M, clu, approaches, blocks, isTwoMode = NULL, isSym = NULL,
    diag = 1, IM = NULL, EM = NULL, Earr = NULL, justChange = FALSE, 
    rowCluChange = c(0, 0), colCluChange = c(0, 0), sameIM = FALSE, 
    regFun = "max", homFun = "ss", usePreSpecM = NULL, preSpecM = NULL, 
    save.initial.param = TRUE, relWeights = 1, posWeights = 1, 
    blockTypeWeights = 1, combWeights = NULL, returnEnv = FALSE)
    
    
optParC(M, nMode = NULL,isSym = NULL, diag = 1, clu, approaches, blocks,
    useMulti = FALSE, maxPar = 50, IM = NULL, EM = NULL, Earr = NULL,
    justChange = FALSE, sameIM = FALSE, regFun = "max", homFun = "ss",
    usePreSpecM = NULL, preSpecM = NULL, minUnitsRowCluster = 1,
    minUnitsColCluster = 1, maxUnitsRowCluster = 9999,
    maxUnitsColCluster = 9999, relWeights = 1, posWeights = 1,
    blockTypeWeights = 1, combWeights = NULL, exchageClusters = "all",
    save.initial.param = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation. The network can have one or more modes (diferent kinds of units with no ties among themselves). If the network is not two-mode, the matrix must be square.}
\item{clu}{A partition. Each unique value represents one cluster. If the nework is one-mode, than this should be a vector, else a list of vectors, one for each mode.}
\item{approaches}{One of the approaches (for each relation in multi-relational netowrks in a vector) described in Žiberna (2006). Possible values are:\cr
  "bin" - binary blockmodeling,\cr
  "val" - valued blockmodeling,\cr
  "ss" - sum of squares homogenity blockmodeling, and\cr
  "ad" - absolute deviations homogenity blockmodeling.}
\item{blocks}{A vector, a list of vectors or an array with names of allowed blocy types. \cr
\cr
Only listing of allowed block types (blockmodel is not pre-specified).\cr
A vector with names of allowed blocktypes. For multi-relational networks, it can be a list of such vectors. For \code{approaches = "bin"} or \code{approaches = "val"}, at least two should be selected. Possible values are:\cr
  \code{"nul"} - null or empty block\cr
  \code{"com"} - complete block\cr
  \code{"rdo"}, \code{"cdo"} - row and column-dominant blocks (binary and valued approach only)\cr
  \code{"reg"} - (f-)regular block\cr
  \code{"rre"}, \code{"cre"} - row and column-(f-)regular blocks\cr
  \code{"rfn"}, \code{"cfn"} - row and column-dominant blocks (binary, valued only)\cr
  \code{"den"} - density block (binary approach only)\cr
  \code{"avg"} - average block (valued approach only)\cr
  \code{"dnc"} - do not care block - the error is always zero\cr
  The ordering is important, since if several block types have identical error, the first on the list is selected.\cr\cr
 A pre-specified blockmodel.\cr
 An array with dimensions four dimensions (see example below). The third and the fourth represent the clusters (for rows and columns). The first is as long as the maximum number of allows block types for a given block. If some block has less possible block types, the empty slots should have values \code{NA}. The second dimension is the number of relations (1 for single-relational networks). The values in the array should be the ones from above.}
\item{isTwoMode}{\code{1} for one-mode networks and \code{2} for two-mode networks. The default value is set to \code{NULL}.}
\item{isSym}{Specifying if the matrix (for each relation) is symetric.}
\item{diag}{Should the special stauts of diagonal be acknowladged. The default value is set to \code{1}.}
\item{IM}{The obtained image for objects.}
\item{EM}{Block errors by blocks.}
\item{Earr}{The array of errors for all allowed block types by next dimensions: allowed block types, relations, row clusters and column clusters. The dimensions should match the dimensions of the block argument if specified as an array.}
\item{justChange}{Value specifying if only the errors for changed clusters should be computed. Used only for debugging purposes by developers.}
\item{rowCluChange}{An array holding the two row clusters where the change occured. Used only for debugging purposes by developers.}
\item{colCluChange}{An array holding the col row clusters where the change occured. Used only for debugging purposes by developers.}
\item{sameIM}{Should we damand the same blockmodel image for all relations. The default value is set to \code{FALSE}.}
\item{regFun}{Function f used in row-f-regular, column-f-regular, and f-regular blocks. Not used in binary approach. For multi-relational networks, it can be a vector of such character strings. The default value is set to \code{"max"}.}
\item{homFun}{In case of homogenity blockmodeling two vairability criteria can be used: \code{"ss"} - sum of squares (set by default) and \code{"ad"} -
absolute deviations.}
\item{usePreSpecM}{Specifiying weather a pre-specified value should be used when computing inconsistency.}
 \item{preSpecM}{Suficient value for individual cells for valued approach. Can be a number or a character string giving the name of a function. Set to \code{"max"} for implicit approach. For multi-relational networks, it can be a vector of such values. In case ob binary blockmodeling this argument is a threshold used for binerizing the network. Therefore all values with values lower than \code{preSpecM} are recoded into 0s, all other into 1s. For multi-relational networks, it can be a vector of such values.}
\item{save.initial.param}{Should the inital parameters (\code{approaches}, ...) be saved. The default value is \code{TRUE}.}
\item{relWeights}{Weights for all type of relations in a blockmodel. The default value is set to \code{1}.}
\item{posWeights}{Weigths for positions in the blockmodel (the dimensions must be the same as the error matrix). For now this is a matix (two-dimensional) even for multi-relational networks.}
\item{blockTypeWeights}{Weights for each type of block used, if they are to be different accros block types (see \code{blocks} above). It must be suplied in form of a named vetor. The default value is set to \code{1}.}
\item{combWeights}{Weights for all type of block used, The default value is set to \code{NULL}.}
\item{returnEnv}{Should the function also return the environment after its completion.}
\item{useMulti}{Which version of local search should be used. The default value is set to \code{FALSE}. If \code{FALSE}, first possible all moves in random order and then all possible exchanges in random order are tired. When a move with lower value of criterion function is found, the algorithm moves to this new partition. If TRUE the version of local search where all possible moves and exchanges are tired first and then the one with the lowest error is selected and used. In this case, several optimal partitions are found. \code{maxPar} best partitions are returned.}
\item{maxPar}{The number of partitions with optimal criterion fuction to be returned. Only used If \code{useMulti} is \code{TRUE}.}
\item{nMode}{Number of nodes. If \code{NULL}, then determined from \code{clu}.}
\item{minUnitsRowCluster}{Minimum number of units in row cluster.}
\item{minUnitsColCluster}{Minimum number of units in col cluster.}
\item{maxUnitsRowCluster}{Maximum number of units in row cluster.}
\item{maxUnitsColCluster}{Maximum number of units in col cluster.}
\item{exchageClusters}{A matrix of dimensions "number of clusters" x "number of clusters" indicating to which clusters can units from a specific cluster be moved. Useful for multilevel blockmodeling or/in some other cases where some units cannot mix.}

  }

%\details{
%}
\value{
  A list:
  \item{M}{The matrix of the network analyzed.}
  \item{err}{The error or inconsistency emplirical network with the ideal network for a given blockmodel (model,approach,...) and paritition.}
  \item{clu}{The analyzed partition.}
  \item{EM}{Block errors by blocks.}
  \item{IM}{The obtained image for objects.}
  \item{BM}{Block means by block - only for Homogeneity blockmodeling.}
  \item{Earr}{The array of errors for all allowed block types by next dimensions: allowed block types, relations, row clusters and column clusters. The dimensions should match the dimensions of the block argument if specified as an array.}
}
\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2014): Blockmodeling of multilevel networks. Social Networks, Oct. 2014, vol. 39, no. 1, 46-61. \url{https://doi.org/10.1016/j.socnet.2014.04.002}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, \enc{Anuška}{Anuska} (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.
}

\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{optRandomParC}}, \code{\link{IM}}, \code{\link{clu}}, \code{\link{err}}, \code{\link{plot.crit.fun}}}
\examples{
## Generating a simple network corresponding to the simple Sum of squares
## Structural equivalence with blockmodel:
## nul com
## nul nul
n <- 20
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(5, 15))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

## Computation of criterion function with the correct partition
res <- critFunC(M = net, clu = clu, approaches = "hom", homFun = "ss", blocks = "com")
res$err ## the error is relativly small
res$BM  ## The block means are around 0 or 4
plot(res)

## Computation of criterion function with the correct partition and correct pre-specified blockmodel
## Prespecified blockmodel used
## nul   com
## nul   nul
B <- array(NA, dim = c(1, 1, 2, 2))
B[1, 1, , ] <- "nul"
B[1, 1, 1, 2] <- "com"
B[1, 1, , ]
res <- critFunC(M = net, clu = clu, approaches = "hom", homFun = "ss", blocks = B)
res$err ## the error is relativly small
res$IM 
plot(res)

## Computation of criterion function with the correct partition 
# and pre-specified blockmodel with some alternatives
## Prespecified blockmodel used
## nul nul|com
## nul nul
B <- array(NA, dim = c(2, 2, 2))
B[1, , ] <- "nul"
B[2, 1, 2] <- "com"
res <- critFunC(M = net, clu = clu, approaches = "hom", homFun = "ss", blocks = B)
res$err ## the error is relativly small
res$IM 
plot(res)

## Computation of criterion function with random partition
clu.rnd <- sample(1:2, size = n, replace = TRUE)
res.rnd <- critFunC(M = net, clu = clu.rnd, approaches = "hom", 
homFun = "ss", blocks = "com")
res.rnd$err ## the error is larger
res.rnd$BM  ## random block means
plot(res.rnd)

#adapt network for Valued blockmodeling with the same model
net[net > 4] <- 4
net[net < 0] <- 0

## Computation of criterion function with the correct partition
res <- critFunC(M = net, clu = clu, approaches = "val", 
blocks = c("nul", "com"), preSpecM = 4)
res$err ## the error is relativly small
res$IM 
## The image corresponds to the one used for generation of
## The network
plot(res)

## Computation of criterion function with random partition
res.rnd <- critFunC(M = net, clu = clu.rnd, approaches = "val",
blocks = c("nul", "com"), preSpecM = 4)
res.rnd$err ## the error is larger
res.rnd$IM ## all blocks are probably nul
plot(res.rnd)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
