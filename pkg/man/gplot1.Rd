\name{gplot1}
\alias{gplot1}
\alias{gplot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A wrapper for funcction gplot - Two-Dimensional Visualization of Graphs}
\description{
The function calls function \code{gplot} from library sna with different defaults. Usefun for ploting image graphs.
}
\usage{
gplot1(M, diag = TRUE, displaylabels = TRUE, boxed.labels = FALSE,
   loop.cex = 4, arrowhead.cex = NULL, arrowheads.fun = "sqrt",
   edge.lwd = 1, edge.col = "default", rel.thresh = 0.05, ...)
gplot2(M, uselen = TRUE, usecurve = TRUE, edge.len = 0.001,
   diag = TRUE, displaylabels = TRUE, boxed.labels = FALSE,
   loop.cex = 4, arrowhead.cex = 2.5, edge.lwd = 1,
   edge.col = "default", rel.thresh = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matix (array) of a graph or set thereof.  This data may be valued. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{rel.thresh}{ real number indicating the lower relative (compared to the highest value) threshold for tie values.  Only ties of value >\code{thresh} are displayed.  By default, \code{thresh}=0.}
  \item{displaylabels}{ boolean; should vertex labels be displayed? }
  \item{boxed.labels}{ boolean; place vertex labels within boxes? }
  \item{arrowhead.cex}{An expansion factor for edge arrowheads.}
  \item{arrowheads.fun}{A function for scaling arrowheands.}
  \item{loop.cex}{ expansion factor for loops; may be given as a vector, if loops are to be of different sizes. }
  \item{edge.col}{ color for edges; may be given as a vector or adjacency matrix, if edges are to be of different colors. }
  \item{edge.lwd}{ line width scale for edges; if set greater than 0, edge widths are scaled by \code{edge.lwd*dat}.  May be given as a vector or adjacency matrix, if edges are to have different line widths. }
  \item{edge.len}{ if \code{uselen==TRUE}, curved edge lengths are scaled by \code{edge.len}. }
  \item{uselen}{ boolean; should we use \code{edge.len} to rescale edge lengths? }
  \item{usecurve}{ boolean; should we use \code{edge.curve}? }
  \item{\dots}{ additional arguments to \code{\link{plot}} or \code{gplot} from package \code{sna}:\cr\cr
  \bold{\code{mode}}:  the vertex placement algorithm; this must correspond to a \code{gplot.layout} function from package \code{sna}.
  }
}

\value{
 Plots a graph
}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{sna:gplot}}
\keyword{graphs}% at least one, from doc/KEYWORDS
