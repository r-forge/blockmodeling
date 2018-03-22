
\name{plot.mat}
\alias{plot.mat}
\alias{plotMat}
%\alias{plot.mat.nm}
\alias{plotMatNm}
\alias{plot.array}
\alias{plotArray}
\alias{plot.crit.fun}
\alias{plot.opt.par}
\alias{plot.opt.par.mode}
\alias{plot.opt.more.par}
\alias{plot.opt.more.par.mode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for plotting a partitioned matrix (representing the network)}
\description{
The main function \code{plot.mat} or \code{plotMat} plots a (optionally partitioned) matrix. If the matrix is partitioned, the rows and columns of the matrix are rearranged according to the partitions. Other functions are only wrappers for \code{plot.mat} or \code{plotMat} for convenience when plotting the results of the corresponding functions. The \code{plotMatNm} plots two matrices based on M, normalized by rows and columns, next to each other. The \code{plot.array} or \code{plotArray} plots an array. \code{plot.mat.nm} has been replaced by \code{plotMatNm}.
}
\usage{
%plot.mat(x = M, M = x, clu = NULL, 
%                ylab = "", xlab = "", main =
%                 NULL, print.val = !length(table(M)) <= 2, print.0 =
%                 FALSE, plot.legend = !print.val && !length(table(M))
%                 <= 2, print.legend.val = "out", print.digits.legend =
%                 2, print.digits.cells = 2, print.cells.mf = NULL,
%                 outer.title = FALSE, title.line = ifelse(outer.title,
%                 -1.5, 7), mar = c(0.5, 7, 8.5, 0) + 0.1, cex.val =
%                 "default", val.y.coor.cor = 0, val.x.coor.cor = 0,
%                 cex.legend = 1, legend.title = "Legend", cex.axes =
%                 "default", print.axes.val = NULL, print.x.axis.val =
%                 !is.null(colnames(M)), print.y.axis.val =
%                 !is.null(rownames(M)), x.axis.val.pos = 1.01,
%                 y.axis.val.pos = -0.01, cex.main = par()$cex.main,
%                 cex.lab = par()$cex.lab, yaxis.line = -1.5, xaxis.line
%                 = -1, legend.left = 0.4, legend.up = 0.03, legend.size
%                 = 1/min(dim(M)), legend.text.hor.pos = 0.5,
%                 par.line.width = 3, par.line.col = "blue", IM.dens =
%                 NULL, IM = NULL, wnet = NULL, wIM = NULL, use.IM =
%                 length(dim(IM)) == length(dim(M)) | !is.null(wIM),
%                 dens.leg = c(null = 100, nul = 100), blackdens = 70,
%                 plotLines = FALSE, frameMatrix = TRUE, x0ParLine =
%                 -0.1, x1ParLine = 1, y0ParLine = 0, y1ParLine = 1.1,
%                 colByUnits = NULL, colByRow = NULL, colByCol = NULL,
%                 mulCol = 2, joinColOperator = "+", colTies = FALSE,
%                 maxValPlot = NULL, printMultipliedMessage = TRUE,
%                 replaceNAdiagWith0 = TRUE, ...)
                
                 
plotMat(x = M, M = x, clu = NULL, 
                ylab = "", xlab = "", main =
                 NULL, print.val = !length(table(M)) <= 2, print.0 =
                 FALSE, plot.legend = !print.val && !length(table(M))
                 <= 2, print.legend.val = "out", print.digits.legend =
                 2, print.digits.cells = 2, print.cells.mf = NULL,
                 outer.title = FALSE, title.line = ifelse(outer.title,
                 -1.5, 7), mar = c(0.5, 7, 8.5, 0) + 0.1, cex.val =
                 "default", val.y.coor.cor = 0, val.x.coor.cor = 0,
                 cex.legend = 1, legend.title = "Legend", cex.axes =
                 "default", print.axes.val = NULL, print.x.axis.val =
                 !is.null(colnames(M)), print.y.axis.val =
                 !is.null(rownames(M)), x.axis.val.pos = 1.01,
                 y.axis.val.pos = -0.01, cex.main = par()$cex.main,
                 cex.lab = par()$cex.lab, yaxis.line = -1.5, xaxis.line
                 = -1, legend.left = 0.4, legend.up = 0.03, legend.size
                 = 1/min(dim(M)), legend.text.hor.pos = 0.5,
                 par.line.width = 3, par.line.col = "blue", IM.dens =
                 NULL, IM = NULL, wnet = NULL, wIM = NULL, use.IM =
                 length(dim(IM)) == length(dim(M)) | !is.null(wIM),
                 dens.leg = c(null = 100, nul = 100), blackdens = 70,
                 plotLines = FALSE, frameMatrix = TRUE, x0ParLine =
                 -0.1, x1ParLine = 1, y0ParLine = 0, y1ParLine = 1.1,
                 colByUnits = NULL, colByRow = NULL, colByCol = NULL,
                 mulCol = 2, joinColOperator = "+", colTies = FALSE,
                 maxValPlot = NULL, printMultipliedMessage = TRUE,
                 replaceNAdiagWith0 = TRUE, ...)


% plot.mat.nm(x = M, M = x, ..., main.title = NULL,
%   title.row = "Row normalized",
%   title.col = "Column normalized",
%   main.title.line = -2, par.set = list(mfrow = c(1, 2)))


plotMatNm(x = M, M = x, ..., main.title = NULL,
   title.row = "Row normalized",
   title.col = "Column normalized",
   main.title.line = -2, par.set = list(mfrow = c(1, 2)))



%plot.array(x = M, M = x, ..., main.title = NULL, main.title.line
%                 = -2, mfrow = NULL)

        
plotArray(x = M, M = x, ..., main.title = NULL, main.title.line
                 = -2, mfrow = NULL)

\method{plot}{mat}(x = M, M = x, clu = NULL, 
                  ylab = "", xlab = "", main =
                 NULL, print.val = !length(table(M)) <= 2, print.0 =
                 FALSE, plot.legend = !print.val && !length(table(M))
                 <= 2, print.legend.val = "out", print.digits.legend =
                 2, print.digits.cells = 2, print.cells.mf = NULL,
                 outer.title = FALSE, title.line = ifelse(outer.title,
                 -1.5, 7), mar = c(0.5, 7, 8.5, 0) + 0.1, cex.val =
                 "default", val.y.coor.cor = 0, val.x.coor.cor = 0,
                 cex.legend = 1, legend.title = "Legend", cex.axes =
                 "default", print.axes.val = NULL, print.x.axis.val =
                 !is.null(colnames(M)), print.y.axis.val =
                 !is.null(rownames(M)), x.axis.val.pos = 1.01,
                 y.axis.val.pos = -0.01, cex.main = par()$cex.main,
                 cex.lab = par()$cex.lab, yaxis.line = -1.5, xaxis.line
                 = -1, legend.left = 0.4, legend.up = 0.03, legend.size
                 = 1/min(dim(M)), legend.text.hor.pos = 0.5,
                 par.line.width = 3, par.line.col = "blue", IM.dens =
                 NULL, IM = NULL, wnet = NULL, wIM = NULL, use.IM =
                 length(dim(IM)) == length(dim(M)) | !is.null(wIM),
                 dens.leg = c(null = 100, nul = 100), blackdens = 70,
                 plotLines = FALSE, frameMatrix = TRUE, x0ParLine =
                 -0.1, x1ParLine = 1, y0ParLine = 0, y1ParLine = 1.1,
                 colByUnits = NULL, colByRow = NULL, colByCol = NULL,
                 mulCol = 2, joinColOperator = "+", colTies = FALSE,
                 maxValPlot = NULL, printMultipliedMessage = TRUE,
                 replaceNAdiagWith0 = TRUE, ...)

\method{plot}{crit.fun}(x, main = NULL, ...)


\method{plot}{array}(x = M, M = x, ..., main.title = NULL, main.title.line
                 = -2, mfrow = NULL)

\method{plot}{opt.par}(x, main = NULL, which = 1, ...)

\method{plot}{opt.par.mode}(x, main = NULL, which = 1, ...)

\method{plot}{opt.more.par}(x, main = NULL, which = 1, ...)

\method{plot}{opt.more.par.mode}(x, main = NULL, which = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A result from a corespodning function or a matrix or similar object representing a network.}
  \item{M}{A matrix or similar object representing a network - either \code{x} or \code{M} must be supplied - both are here to make the code compatible with generic and with older functions.}
  \item{clu}{A partition. Each unique value represents one cluster. If the network is one-mode,
than this should be a vector, else a list of vectors, one for each mode.}
  \item{ylab}{Label for y axis.}
  \item{xlab}{Label for x axis.}
  \item{main}{Main title.}
  \item{main.title}{Main title in plot.array version.}
  \item{main.title.line}{The line in which main title is printed in \code{plot.array} version.}
  \item{mfrow}{\code{mfrow} argument to \code{par} - number of row and column plots to be ploted on one figure.}
  \item{print.val}{Should the values be printed in the matrix.}
  \item{print.0}{If \code{print.val = TRUE} Should the 0s be printed in the matrix.}
  \item{plot.legend}{Should the legend for shades be ploted.}
  \item{print.legend.val}{Should the values be printed in the legend.}
  \item{print.digits.legend}{The number of digits that should appear in the legend.}
  \item{print.digits.cells}{The number of digits that should appear in the cells (of the matrix and/or legend).}
  \item{print.cells.mf}{if not \code{NULL}, the above argument is igonred, the cell values are printed as the cell are multiplied by this factor and rounded.}
  \item{outer.title}{Should the title be printed on the 'inner' or 'outer' plot, default is \code{'inner'}.}
  \item{title.line}{The line (from the top) where the title should be printed. The suitable values depend heavily on the displey type.}
  \item{mar}{A numerical vector of the form \code{c(bottom, left, top, right)} which gives the lines of margin to be specified on the four sides of the plot. The R default for ordianry plots is \code{c(5, 4, 4, 2) + 0.1}, while this functions default is \code{c(0.5, 7, 8.5, 0) + 0.1}.}
  \item{cex.val}{Size of the values printed. The default is \code{10/'number of units'}.}
  \item{val.y.coor.cor}{Correction for centering the values in the sqares in y direction.}
  \item{val.x.coor.cor}{Correction for centering the values in the sqares in x direction.}
  \item{cex.legend}{Size of the text in the legend.}
  \item{legend.title}{The title of the legend.}
  \item{cex.axes}{Size of the characters in axes. Default makes the cex so small that all categories can be printed.}
  \item{print.axes.val}{Should the axes values be printed. Default prints each axis if \code{rownames} or \code{colnames} is not \code{NULL}.}
  \item{print.x.axis.val}{Should the x axis values be printed. Default prints each axis if \code{rownames} or \code{colnames} is not \code{NULL}.}
  \item{print.y.axis.val}{Should the y axis values be printed. Default prints each axis if \code{rownames} or \code{colnames} is not \code{NULL}.}
  \item{x.axis.val.pos}{x coordiante of the y axis values.}
  \item{y.axis.val.pos}{y coordiante of the x axis values.}
  \item{cex.main}{Size of the text in the main title.}
  \item{cex.lab}{Size of the text in matrix.}
  \item{yaxis.line}{The position of the y axis (the argument 'line').}
  \item{xaxis.line}{The position of the x axis (the argument 'line').}
  \item{legend.left}{How much left should the legend be from the matrix.}
  \item{legend.up}{How much up should the legend be from the matrix.}
  \item{legend.size}{Relative legend size.}
  \item{legend.text.hor.pos}{Horizontal position of the legend text (bottom) - 0 = bottom, 0.5 = middle,...}
  \item{par.line.width}{The width of the line that seperates the partitions.}
  \item{par.line.col}{The color of the line that seperates the partitions.}
  \item{IM.dens}{The densitiey of shading lines for each block.}
  \item{IM}{The image (as obtaind with \code{critFunC}) of the blockmodel. \code{dens.leg} is used to translate this image into \code{IM.dens}.}
   \item{wnet}{Specifies which matrix (if more) should be ploted - used if \code{M} is an array.}
   \item{wIM}{Specifies which IM (if more) should be used for ploting. The default value si set to \code{wnet}) - used if IM is an array.}
  \item{use.IM}{Specifies if IM should be used for plotting.}
    \item{dens.leg}{It is used to translate the \code{IM} into \code{IM.dens}.}
  \item{blackdens}{At which density should the values on dark colurs of lines be printed in white.}
  \item{plotLines}{Should the lines in the matrix be printed. The default value is set to \code{FALSE}, best set to \code{TRUE} for very small networks.}
  \item{frameMatrix}{Should the matrix be framed (if plotLines is \code{FALSE}). The default value is set to \code{TRUE}.}
  \item{x0ParLine}{Coordinates for lines seperating clusters.}
  \item{x1ParLine}{Coordinates for lines seperating clusters.}
  \item{y0ParLine}{Coordinates for lines seperating clusters.}
  \item{y1ParLine}{Coordinates for lines seperating clusters.}
   \item{colByUnits}{Coloring units units. It should be a vector of unit length.}
  \item{colByRow}{Coloring units by rows. It should be a vector of unit length.}
  \item{colByCol}{Coloring units by columns. It should be a vector of unit length.}
  \item{mulCol}{Multiply color when joining with row, column. Only used when when \code{colByUnits} is not NULL.}
  \item{joinColOperator}{Function to join \code{colByRow} and \code{colByCol}. The default value is set to \code{"+"}.}
  \item{colTies}{If \code{TRUE}, ties are colored, if \code{FALSE}, 0-ties are colored.}
  \item{maxValPlot}{The value to use as maximum when computing colors (ties with maximal positive value are ploted as black).}
  \item{printMultipliedMessage}{Should the message '* all values in cells were multiplied by' be printed on the plot. The default value is set to \code{TRUE}.}
  \item{replaceNAdiagWith0}{If \code{replaceNAdiagWith0 = TRUE} Should the NA values on diagonal of a matrix be replaced with 0s}
   \item{title.row}{Title for the row-normalized matrix in nm version}
   \item{title.col}{Title for the column-normalized matrix in nm version}
   \item{par.set}{A list of possible ploting paramters (to \code{par}) to be used in nm version}
  \item{which}{Which (if there are more than one) of optimal solutions to plot.}
  \item{\dots}{Aditional arguments to \code{plot.default} for \code{plotMat} and also to \code{plotMat} for other functions.}
}
\value{
  The functions are used for their side affect - plotting.
}
\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{critFunC}}, \code{\link{optRandomParC}}}
\examples{
## Generation of the network
n <- 20
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(5, 15))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

## Ploting the network
plotMat(M = net, clu = clu, print.digits.cells = 3)
class(net) <- "mat"
plot(net, clu = clu)
## See corespodning functions for examples for other plotting
## functions
## presented, that are esentially only the wrappers for "plot.max"
}
\keyword{graphs}% at least one, from doc/KEYWORDS
\keyword{hplot}% at least one, from doc/KEYWORDS
