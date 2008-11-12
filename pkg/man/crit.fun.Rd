\name{crit.fun}
\alias{crit.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the criterion function for a given network and partition}
\description{
The function computes the value of a criterion function for a given network and partition for Generalized blockmodeling. (Žiberna, 2006) based on other parameters (see below).
}
\usage{
crit.fun(M, clu, approach, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{clu}{A partition. Each unique value represents one cluster. If the nework is one-mode, than this should be a vector, else a list of vectors, one for each mode}
  \item{approach}{One of the approaches (for each relation in multi-relational netowrks in a vector) described in Žiberna (2006). Possible values are:\cr
  "bin" - binary blockmodeling,\cr
  "val" - valued blockmodeling,\cr
  "imp" - implicit blockmodeling,\cr
  "ss" - sum of squares homogenity blockmodeling, and\cr
  "ad" - absolute deviations homogenity blockmodeling.}
  \item{\dots}{Several arguments, which are explaind below. They are actually used by the function \code{gen.crit.fun}, however since this functions is not intented to be called directly, it also has no help files. Therefore these arguments are described below. Which are needed depends on the \code{approach} selected}
  \item{blocks}{A vector with names of allowed blocktypes. For multi-relational networks, it can be a list of such vectors. For approaches "bin", and "val", at least two should be selected. Possible values are are:\cr
  "null" - null or empty block\cr
  "com" - complete block\cr
  "rdo", "cdo" - row and column-dominant blocks (binary, valued, and implicit approach only)\cr
  "reg" - (f-)regular block\cr
  "rre", "cre" - row and column-(f-)regular blocks\cr
  "rfn", "cfn" - row and column-dominant blocks (binary, valued, and implicit approach only)\cr
  "den" - density block (binary approach only)\cr
  "avg" - average block (valued and implicit approach only)\cr
  "dnc" - do not care block - the error is always zero\cr
  The ordering is important, since if several block types have identical error, the first on the list is selected.}
  \item{BLOCKS}{An alternative to \code{blocks}. A pre-specified blockmodel. An array with dimensions three dimensions (see example below). The second and the third repreent the clusters (for rows and columns), while the first is as long as the maxsimum number of allows block types for a given block. If some block has less possible block types, the empty slots should have values \code{NA}. The values in the array should be the ones from above. For multi-relational networks, it can be a list of such arrays.}
  \item{m}{Suficient value for individual cells for valued approach. Can be a number or a character string giving the name of a function. Set to \code{"max"} for implicit approach. For multi-relational networks, it can be a vector of such values.}
  \item{cut}{(default = \code{min(M[M > 0])}) The threshold used for binerizing the network for use with binary blockmodeling. All values with values lower than \code{cut} are recoded into 0s, all other into 1s. For multi-relational networks, it can be a vector of such values.}
  \item{FUN}{(default = "max") Function f used in row-f-regular, column-f-regular, and f-regular blocks. Not used in binary approach. For multi-relational networks, it can be a vector of such character strings.}
  \item{norm}{Should the block errors (inconsistencies) be normalized with the size of the blocks, the block error does not depend on block size? The default is \code{FALSE}. Original version of implicit approach suggests \code{TRUE}, however the default is \code{FALSE} even for this approach based on better results in simulations. For multi-relational networks, it can be a vector.}
  \item{normMto2rel} {Create two-realation netowrk from one relational network through row and column normalization.  The default is \code{FALSE}}
  \item{sameModel} {Should we damand the same blockmodel for all relations. If set to TRUE, it demands that accros all relations the ideal block on the same position in the matrix BLOCKS should be chosen. Usually, these positions are occupied by the same blocks. If not, use with caution. The default is the value of \code{normMto2rel}.}
  \item{normbym}{ The default is \code{FALSE} for valued and implicit approach, elsewher not used. Original version of implicit approach suggests \code{TRUE}, however the default is \code{FALSE} even for this approach based on better results in simulations. For multi-relational networks, it can be a vector.}
  \item{allow.max0}{Should the maximum that is the basis for calculation of inconsistencies in implicit blockmodeling be allowed to be 0. If FALSE, the maximum is in such case set to the maximum of the network (if maximum of a block is 0) or to the maxsimum of the block (if row or column maxismum is 0)Used only in implicit blockmodeling. If \code{TRUE}, the incosistency of an ideal null block is 0 for all block types. The default is \code{FALSE} if null blocks are incluede in the allowed blocks in at least one block and \code{TRUE} otherwise.}
  \item{allow.dom0}{Should the dominant row or column (in row- or column-dominant blocks) be allowed to be 0. Used only in implicit blockmodeling. The default is \code{FALSE}.}
  \item{max.con.val} {Should the largest values be cencored, limited to (larger values set to) - resonoble values are:\cr
  "non" - (the default) no transformation is done\cr
  "m" - (the default for implicit blockmodeling) the maximum value equals the value of the parameter m\cr
  numerical values (usually) larger then parameter m and lower the the maximum value in M}
  \item{mindim}{(default = 2) Minimal dimension (number of rows or columns) demanded for row and column-dominant and -functional blocks.}
  \item{mindimreg}{(default = \code{FALSE}) Should the mindim argument also be used for (row or coulum-)(f-)regular blocks}
  \item{blockWeights}{Weights for each type of block used, if they are to be different accros block types (see \code{blocks} above). It must be suplied in form of a named vetor \code{blockWeights = c(name.of.block.type1=weight,...)}. If some of the block types used are not listed, they are given weight 1.}
  \item{positionWeights}{weigths for positions in the blockmodel (the dimensions must be the same as the error matrix). For now this is a matix (two-dimensional) even for multi-relational networks.}
  \item{save.err.v}{(default = \code{FALSE}) Should the error vector for all allowed block types in each block be saved?}
  \item{BLOCK.CV}{An array with the same dimmesions as  \code{BLOCKS} of central values  for pre-specified homogenity (sum of squares and absolute deviations) approach. For multi-relational networks, it can be a list of such arrays.}
  \item{CV.use}{An array with the same dimmesions as  \code{BLOCKS.CV} with instuctions how to treat these centarl values. For multi-relational networks, it can be a list of such arrays. Possible alternatives are:\cr
  "fixed" - the central value is fixed to the value specified in \code{BLOCKS.CV}.\cr
  "min" - the central value specified in \code{BLOCKS.CV} is the minimal possible central value for a block. The central value for the block is computed as the maximum of the value specified in \code{BLOCKS.CV} and the empirical value computed based on tie values in the block.\cr
  "max" - the central value specified in \code{BLOCKS.CV} is the maximal possible central value for a block. The central value for the block is computed as the minimum of the value specified in \code{BLOCKS.CV} and the empirical value computed based on tie values in the block.\cr
  "free" - the central value is free, the value specified in \code{BLOCKS.CV} is igneored. The central value for the block is computed as the empirical value computed based on tie values in the block.\cr}
  \item{use.for}{(default = \code{TRUE}) Should FORTRAN subrutines be used where available (available for only very special cases, currently only for using "ss" aproach and only complete blocks. If you are using such setting and some special features (these are not implemented in FORTRAN subrutines), it's safer to set it to FASLE, as the fuction may miss that these features are not implemented in FORTRAN subrutines and use them nevertheless, leading to wrong results.}
  \item{diag}{(default = \code{TRUE}) Should the special stauts of diagonal be acknowladged.}
}
%\details{
%}
\value{
  A list:
  \item{M}{The matrix of the network analyzed}
  \item{err}{The error or inconsistency emplirical network with the ideal network for a given blockmodel (model,approach,...) and paritition}
  \item{clu}{The analyzed partition}
  \item{E}{Block errors by blocks}
  \item{IM}{The obtained image}
  \item{BM}{Block means by block - only for Homogeneity blockmodeling}
  \item{ERR.V}{If selected. The error vector of errors for all allowed block types by blocks. The dimmensions are [rows, columns (,relations - if more than 1)]. Each cell contains a list of errors by block types}
}
\references{
ŽIBERNA, Aleš (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

ŽIBERNA, Aleš. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, Anuška (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.
}

\author{Aleš Žiberna}
\seealso{\code{\link{opt.par}},\code{\link{opt.random.par}},\code{\link{opt.these.par}},\code{\link{check.these.par}},\code{\link{plot.crit.fun}}}
\examples{
#generating a simple network corresponding to the simple Sum of squares structural equivalence with blockmodel:
# null com
# null null
n<-20
net<-matrix(NA,ncol=n,nrow=n)
clu<-rep(1:2,times=c(5,15))
tclu<-table(clu)
net[clu==1,clu==1]<-rnorm(n=tclu[1]*tclu[1],mean=0,sd=1)
net[clu==1,clu==2]<-rnorm(n=tclu[1]*tclu[2],mean=4,sd=1)
net[clu==2,clu==1]<-rnorm(n=tclu[2]*tclu[1],mean=0,sd=1)
net[clu==2,clu==2]<-rnorm(n=tclu[2]*tclu[2],mean=0,sd=1)

#computation of criterion function with the correct partition
res<-crit.fun(M=net,clu=clu,approach="ss",blocks="com")
res$err #the error is relativly small
res$BM  #The block means are around 0 or 4
plot(res)

#computation of criterion function with random partition
clu.rnd<-sample(1:2,size=n,replace=TRUE)
res.rnd<-crit.fun(M=net,clu=clu.rnd,approach="ss",blocks="com")
res.rnd$err #the error is larger
res.rnd$BM  #random block means
plot(res.rnd)

#adapt network for Valued blockmodeling with the same model
net[net>4]<-4
net[net<0]<-0

#computation of criterion function with the correct partition
res<-crit.fun(M=net,clu=clu,approach="val",blocks=c("null","com"),m=4)
res$err #the error is relativly small
res$IM  #The image corresponds to the one used for generation of the network
plot(res)

#computation of criterion function with random partition
res.rnd<-crit.fun(M=net,clu=clu.rnd,approach="val",blocks=c("null","com"),m=4)
res.rnd$err #the error is larger
res.rnd$IM #all blocks are probably null
plot(res.rnd)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{graphs}% at least one, from doc/KEYWORDS
