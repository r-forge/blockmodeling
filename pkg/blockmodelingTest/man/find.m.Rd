\name{find.m}
\alias{find.m}
\alias{find.m2}
\alias{find.cut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computing the threshold}
\description{
The functions compute the maxsimum value of \code{m/cut} where a ceratin block is still classified as \code{alt.blocks} and not "null".
The difference between \code{find.m} and \code{find.m2} it that \code{find.m} uses an optimizational approach and is faster and more precise than \code{find.m2}. However, \code{find.m} only supports regular ("reg") and complete ("com") as \code{alt.blocks}, while \code{find.m2} supports all block types. Also, \code{find.m} does not always work, sepecially if \code{cormat} is not "none".
}
\usage{
find.m(M, clu, alt.blocks = "reg", diag = !is.list(clu),
    cormet = "none", half = TRUE, FUN = "max")
find.m2(M, clu, alt.blocks = "reg", neval = 100, half = TRUE,
    ms = NULL, ...)
find.cut(M, clu, alt.blocks = "reg", cuts = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network can have one or more modes (diferent kinds of units with no ties among themselvs. If the network is not two-mode, the matrix must be square.}
  \item{clu}{A partition. Each unique value represents one cluster. If the nework is  one-mode, than this should be a vector, else a list of vectors, one for each mode.}
  \item{alt.blocks}{Only one of allowed blocktypes, as alternative to the null block:\cr
  "com" - complete block\cr
  "rdo", "cdo" - row and column-dominant blocks (binary, valued, and implicit approach only)\cr
  "reg" - (f-)regular block\cr
  "rre", "cre" - row and column-(f-)regular blocks\cr
  "rfn", "cfn" - row and column-dominant blocks (binary, valued, and implicit approach only)\cr
  "den" - density block (binary approach only)\cr
  "avg" - average block (valued approach only)}
  \item{diag}{(default = \code{TRUE}) Should the special stauts of diagonal be acknowladged.}
  \item{cormet}{Which metho should be used to correct for diferent maxismum error contributins?\cr
      "none" - no correction\cr
      "censor" - censor values larger than m\cr
      "correct" -  so that the maxsimum possible error contribution of the cell is the same regardles of a condition (either that somthing must be o or at least \code{m}).}
  \item{FUN}{(default = "max") Function f used in row-f-regular, column-f-regular, and f-regular blocks.}
  \item{cuts}{The cuts which should be evaluatated. If \code{cuts="all"n} (default), all unique values are evaluated.}
  \item{neval}{Number of different \code{m} values to be evaluated.}
  \item{half}{Should the returned value of m be one half of the value where the incosnistencies are the same.}
  \item{ms}{The values of m where the function should be evaluated.}
  \item{\dots}{Other parameters to \code{crit.fun}.}
}
\value{
  A matrix of maximal \code{m/cut} values.
}
\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.

DOREIAN, Patrick, BATAGELJ, Vladimir, FERLIGOJ, \enc{Anuška}{Anuska} (2005): Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press, 2005. XV, 384 p., ISBN 0-521-84085-6.
}

\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{critFunC}} and maybe also \code{\link{optParC}}, \code{\link{plotMat}}}
\keyword{cluster}% at least one, from doc/KEYWORDS
