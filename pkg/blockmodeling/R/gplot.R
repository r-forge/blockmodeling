#' A wrapper for function  gplot - Two-Dimensional Visualization of Graphs
#' 
#' The function calls function \code{gplot} from the library \code{sna} with different defaults. Use fun for plotting image graphs.
#' 
# #' @usage gplot1(M, diag = TRUE,
# #'     displaylabels = TRUE, boxed.labels = FALSE,
# #'     loop.cex = 4, edge.lwd = 1, edge.col = "default",
# #'     rel.thresh = 0.05, ...)
#' 
# #' gplot2(M, uselen = TRUE, usecurve = TRUE,
# #'     edge.len = 0.001, diag = TRUE, 
# #'     displaylabels = TRUE, boxed.labels = FALSE,
# #'     loop.cex = 4, arrowhead.cex = 2.5, 
# #'     edge.lwd = 1, edge.col = "default", rel.thresh = 0.05, ...)
#'
#' @param M A matrix (array) of a graph or set thereof.  This data may be valued.
#' @param diag Boolean indicating whether or not the diagonal should be treated as valid data
#' Set this \code{TRUE} if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default.
#' @param rel.thresh Real number indicating the lower relative (compared to the highest value) threshold for tie values.
#' Only ties of value \code{thresh} are displayed.  By default, \code{thresh = 0}.
#' @param displaylabels Boolean; should vertex labels be displayed.
#' @param boxed.labels Boolean; place vertex labels within boxes.
#' @param arrowhead.cex An expansion factor for edge arrowheads.
#' @param loop.cex An expansion factor for loops; may be given as a vector, if loops are to be of different sizes.
#' @param edge.lwd Line width scale for edges; if set greater than 0, edge widths are scaled by \code{edge.lwd*dat}.
#' May be given as a vector or adjacency matrix, if edges are to have different line widths.
#' @param edge.col Color for edges; may be given as a vector or adjacency matrix, if edges are to be of different colors.
#' @param edge.len If \code{uselen == TRUE}, curved edge lengths are scaled by \code{edge.len}.
#' @param uselen Boolean; should we use \code{edge.len} to rescale edge lengths.
#' @param usecurve Boolean; should we use \code{edge.curve}.
#' @param \dots Additional arguments to \code{\link{plot}} or \code{link{sna::gplot}}:\cr\cr
#' \bold{\code{mode}}:  the vertex placement algorithm; this must correspond to a \code{gplot.layout} function from package \code{sna}. 
#'
#' @return Plots a graph.
#'
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @seealso \code{link{sna::gplot}}
#' @keywords graphs
#' @importFrom grDevices gray
#' 
#' @export
 
"gplot1" <-function(M,diag=TRUE,displaylabels=TRUE,boxed.labels=FALSE,loop.cex=4,edge.lwd=1,edge.col="default",rel.thresh=0.05,...){
  if(requireNamespace("sna", quietly = TRUE)){
    M[M<(max(M)*rel.thresh)]<-0  
    if(edge.col[1]=="default") edge.col<-gray(1-M/max(M))
    edge.col<-edge.col[edge.col!=gray(1)]
    sna::gplot(dat=M,diag=diag,displaylabels=displaylabels,boxed.labels=boxed.labels,loop.cex=loop.cex,edge.lwd=edge.lwd,edge.col=edge.col,...)
  } else stop("Package \"sna\" is needed for this function to work. Please install it.",
                   call. = FALSE)
}

#' @rdname gplot1
#' 
#' @export

"gplot2" <-
function(M,uselen=TRUE,usecurve=TRUE,edge.len=0.001,diag=TRUE,displaylabels=TRUE,boxed.labels=FALSE,loop.cex=4,arrowhead.cex=2.5,edge.lwd=1,edge.col="default",rel.thresh=0.05,...){
  if(requireNamespace("sna", quietly = TRUE)){
    M[M<(max(M)*rel.thresh)]<-0
    if(edge.col[1]=="default") edge.col<-gray(1-M/max(M))
    edge.col<-edge.col[edge.col!=gray(1)]
    sna::gplot(dat=M,uselen=uselen,usecurve=usecurve,edge.len=edge.len,diag=diag,displaylabels=displaylabels,boxed.labels=boxed.labels,loop.cex=loop.cex,arrowhead.cex=arrowhead.cex,edge.lwd=edge.lwd,edge.col=edge.col,...)
  } else stop("Package \"sna\" is needed for this function to work. Please install it.",
              call. = FALSE)
}

