\name{gplot1}
\alias{gplot1}
\alias{gplot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A wrapper for function  gplot - Two-Dimensional Visualization of Graphs}
\description{
The function calls function \code{gplot} from the library \code{sna} with different defaults. Use fun for plotting image graphs.
}
\usage{
gplot1(M, diag = TRUE, 
    displaylabels = TRUE, boxed.labels = FALSE,
    loop.cex = 4, edge.lwd = 1, edge.col = "default",
    rel.thresh = 0.05, ...)
    
gplot2(M, uselen = TRUE, usecurve = TRUE,
    edge.len = 0.001, diag = TRUE, 
    displaylabels = TRUE, boxed.labels = FALSE,
    loop.cex = 4, arrowhead.cex = 2.5, 
    edge.lwd = 1, edge.col = "default", rel.thresh = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix (array) of a graph or set thereof.  This data may be valued. }
  \item{diag}{Boolean indicating whether or not the diagonal should be treated as valid data.  Set this \code{TRUE} if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default.}
  \item{rel.thresh}{Real number indicating the lower relative (compared to the highest value) threshold for tie values.  Only ties of value \code{thresh} are displayed.  By default, \code{thresh = 0}.}
  \item{displaylabels}{Boolean; should vertex labels be displayed. }
  \item{boxed.labels}{Boolean; place vertex labels within boxes. }
  \item{arrowhead.cex}{An expansion factor for edge arrowheads.}
  \item{loop.cex}{ expansion factor for loops; may be given as a vector, if loops are to be of different sizes. }
  \item{edge.col}{Color for edges; may be given as a vector or adjacency matrix, if edges are to be of different colors. }
  \item{edge.lwd}{Line width scale for edges; if set greater than 0, edge widths are scaled by \code{edge.lwd*dat}.  May be given as a vector or adjacency matrix, if edges are to have different line widths. }
  \item{edge.len}{ if \code{uselen == TRUE}, curved edge lengths are scaled by \code{edge.len}. }
  \item{uselen}{Boolean; should we use \code{edge.len} to rescale edge lengths.}
  \item{usecurve}{Boolean; should we use \code{edge.curve}.}
  \item{\dots}{Additional arguments to \code{\link{plot}} or \code{gplot} from package \code{sna}:\cr\cr
  \bold{\code{mode}}:  the vertex placement algorithm; this must correspond to a \code{gplot.layout} function from package \code{sna}.
  }
}

\value{
 Plots a graph.
}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{sna:gplot}}
\keyword{graphs}% at least one, from doc/KEYWORDS
