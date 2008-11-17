\name{plot.mat}
\alias{plot.mat}
\alias{plot.mat.nm}
\alias{plot.crit.fun}
\alias{plot.opt.par}
\alias{plot.opt.par.mode}
\alias{plot.opt.more.par}
\alias{plot.opt.more.par.mode}
\alias{plot.check.these.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for plotting a partitioned matrix}
\description{
The main function \code{plot.mat} plots a (optionally partitioned) matrix. If the matrix is partitioned, the rows and columns of the matrix are rearranged according to the partitions. Other functions are only wrappers for \code{plot.mat} for convenience when plotting the results of the corresponding functions. The \code{plot.mat.nm} plots two matrices based on M, normalized by rows and columns, next to each other.
}
\usage{
plot.mat(x=M, M=x, clu = NULL, ylab = "", xlab = "", main = NULL,
   print.val = !length(table(M)) <= 2, print.0 = FALSE,
   plot.legend = !print.val && !length(table(M)) <= 2,
   print.legend.val = "out", print.digits.legend = 2,
   print.digits.cells = 2, print.cells.mf = NULL,
   outer.title = !plot.legend,
   title.line = ifelse(outer.title, -1.5, 7),
   mar = c(0.5, 7, 8.5, 0) + 0.1, cex.val = "default",
   val.y.coor.cor = 0, val.x.coor.cor = 0, cex.legend = 1,
   legend.title = "Legend", cex.axes = "default",
   print.axes.val = NULL,
   print.x.axis.val = !is.null(colnames(M)),
   print.y.axis.val = !is.null(rownames(M)),
   x.axis.val.pos = 1.1, y.axis.val.pos = -0.1,
   cex.main = par()$cex.main, cex.lab = par()$cex.lab,
   yaxis.line = -1.5, xaxis.line = -1, legend.left = 0.4,
   legend.up = 0.03, legend.size = 1/min(dim(M)),
   legend.text.hor.pos = 0.5, par.line.width = 3,
   par.line.col = "blue", IM.dens = NULL, IM = NULL, wnet = 1,
   wIM = NULL,
   use.IM = length(dim(IM))==length(dim(M))|!is.null(wIM),
   dens.leg = c(null = 100),
   blackdens = 70, plotLines = TRUE, ...)
  
plot.mat.nm(x=M, M=x, ..., main.title = NULL,
   title.row = "Row normalized", 
   title.col = "Column normalized",
   main.title.line = -2, par.set = list(mfrow = c(1, 2)))

\method{plot}{mat}(x=M, M=x, clu = NULL, ylab = "", xlab = "",
   main = NULL, print.val = !length(table(M)) <= 2,
   print.0 = FALSE,
   plot.legend = !print.val && !length(table(M)) <= 2,
   print.legend.val = "out", print.digits.legend = 2,
   print.digits.cells = 2, print.cells.mf = NULL,
   outer.title = !plot.legend,
   title.line = ifelse(outer.title, -1.5, 7),
   mar = c(0.5, 7, 8.5, 0) + 0.1, cex.val = "default",
   val.y.coor.cor = 0, val.x.coor.cor = 0, cex.legend = 1,
   legend.title = "Legend", cex.axes = "default",
   print.axes.val = NULL,
   print.x.axis.val = !is.null(colnames(M)),
   print.y.axis.val = !is.null(rownames(M)),
   x.axis.val.pos = 1.1, y.axis.val.pos = -0.1,
   cex.main = par()$cex.main, cex.lab = par()$cex.lab,
   yaxis.line = -1.5, xaxis.line = -1, legend.left = 0.4,
   legend.up = 0.03, legend.size = 1/min(dim(M)),
   legend.text.hor.pos = 0.5, par.line.width = 3,
   par.line.col = "blue", IM.dens = NULL, IM = NULL, wnet = 1,
   wIM = NULL,
   use.IM = length(dim(IM)) == length(dim(M)) | !is.null(wIM),
   dens.leg = c(null = 100), blackdens = 70, plotLines = TRUE,
   ...)

\method{plot}{crit.fun}(x, main = NULL, ...)

\method{plot}{opt.par}(x, main = NULL, which = 1, ...)

\method{plot}{opt.par.mode}(x, main = NULL, which = 1, ...)

\method{plot}{opt.more.par}(x, main = NULL, which = 1, ...)

\method{plot}{opt.more.par.mode}(x, main = NULL, which = 1, ...)

\method{plot}{check.these.par}(x, main = NULL, which = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A result from a corespodning function or a matrix or similar object representing a network}
  \item{M}{A matrix or similar object representing a network - either x or M must be supplied - both are here to make the code compatible with generic and with older functions}
  \item{clu}{A partition}
  \item{ylab}{Label for y axis}
  \item{xlab}{Label for x axis}
  \item{main}{Main title}
  \item{main.title}{Main title in nm version}
  \item{main.title.line}{The line in which main title is printed in nm version}
  \item{title.row}{Title for the row-normalized matrix in nm version}
  \item{title.col}{Title for the column-normalized matrix in nm version}
  \item{par.set}{A list of possible ploting paramters (to \code{par}) to be used in nm version}
  \item{print.val}{Should the values be printed in the matrix}
  \item{print.0}{If \code{print.val=TRUE} Should the 0s be printed in the matrix}
  \item{plot.legend}{Should the legend for shades be ploted}
  \item{print.legend.val}{Should the values be printed in the legend}
  \item{print.digits.legend}{The number of digits that should appear in the legend}
  \item{print.digits.cells}{The number of digits that should appear in the cells (of the matrix and/or legend)}
  \item{print.cells.mf}{if not \code{NULL}, the above argument is igonred, the cell values are printed as the cell are multiplied by this factor and rounded}
  \item{outer.title}{Should the title be printed on the 'inner' or 'outer' plot, default is 'inner' if legend is ploted and 'outer' otherwise. May be soon omited.}
  \item{title.line}{The line (from the top) where the title should be printed. The suitable values depend heavily on the displey type.}
  \item{mar}{A numerical vector of the form 'c(bottom, left, top, right)' which gives the lines of margin to be specified on the four sides of the plot. The R default for ordianry plots is 'c(5, 4, 4, 2) + 0.1', while this functions default is c(0.5, 7, 8.5, 0) + 0.1.}
  \item{cex.val}{Size of the values printed. The "default" is 10/"number of units"}
  \item{val.y.coor.cor}{Correction for centering the values in the sqares in y direction}
  \item{val.x.coor.cor}{Correction for centering the values in the sqares in x direction}
  \item{cex.legend}{Size of the text in the legend}
  \item{legend.title}{The title of the legend}
  \item{cex.axes}{Size of the characters in axes, 'default' makes the cex so small that all categories can be printed}
  \item{print.axes.val}{Should the axes values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'}
  \item{print.x.axis.val}{Should the x axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'}
  \item{print.y.axis.val}{Should the y axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'}
  \item{x.axis.val.pos}{x coordiante of the y axis values}
  \item{y.axis.val.pos}{y coordiante of the x axis values}
  \item{cex.main}{Size of the text in the main title}
  \item{cex.lab}{Size of the text in matrix}
  \item{yaxis.line}{The position of the y axis (the argument 'line')}
  \item{xaxis.line}{The position of the x axis (the argument 'line')}
  \item{legend.left}{How much left should the legend be from the matrix}
  \item{legend.up}{How much up should the legend be from the matrix}
  \item{legend.size}{Relative legend size}
  \item{legend.text.hor.pos}{Horizontal position of the legend text (bottom) - 0 = bottom, 0.5 = middle,...}
  \item{par.line.width}{The width of the line that seperates the partitions}
  \item{par.line.col}{The color of the line that seperates the partitions}
  \item{IM.dens}{The densitiey of shading lines for each block}
  \item{IM}{The image (as obtaind with \code{crit.fun}) of the blockmodel. \code{dens.leg} is used to translate this image into \code{IM.dens}.}
  \item{dens.leg}{It is used to translate the \code{IM} into \code{IM.dens}.}
  \item{blackdens}{At which density should the values on dark colurs of lines be printed in white.}
  \item{plotLines}{Should the lines in the matrix be printed - default TRUE, best set to FALSE for larger networks.}
  \item{which}{Which (if there are more than one) of optimal solutions to plot}
  \item{wnet}{Specifies which net (if more) should be ploted - used if M is an array.}
  \item{wIM}{Specifies which IM (if more) should be used for ploting (defualt = wnet) - used if IM is an array.}
  \item{use.IM}{Specifies if IM should IM be used for ploting? be used for ploting?}
  \item{\dots}{Aditional arguments to \code{plot.default} for \code{plot.mat} and also to \code{plot.mat} for other functions}
}
\value{
  The functions are used for their side affect - plotting.
}
\references{
\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales} (2006): Generalized Blockmodeling of Valued Networks. Social Networks, Jan. 2007, vol. 29, no. 1, 105-126. \url{http://dx.doi.org/10.1016/j.socnet.2006.04.002}.

\enc{ŽIBERNA, Aleš}{ZIBERNA, Ales}. Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. J. math. sociol., 2008, vol. 32, no. 1,  57-84. \url{http://www.informaworld.com/smpp/content?content=10.1080/00222500701790207}.
}
\author{\enc{Aleš Žiberna}{Ales Ziberna}}
\seealso{\code{\link{crit.fun}}, \code{\link{opt.par}}, \code{\link{opt.random.par}}, \code{\link{opt.these.par}}, \code{\link{check.these.par}}}
\examples{
#Generation of the network
n<-20
net<-matrix(NA,ncol=n,nrow=n)
clu<-rep(1:2,times=c(5,15))
tclu<-table(clu)
net[clu==1,clu==1]<-rnorm(n=tclu[1]*tclu[1],mean=0,sd=1)
net[clu==1,clu==2]<-rnorm(n=tclu[1]*tclu[2],mean=4,sd=1)
net[clu==2,clu==1]<-rnorm(n=tclu[2]*tclu[1],mean=0,sd=1)
net[clu==2,clu==2]<-rnorm(n=tclu[2]*tclu[2],mean=0,sd=1)

#Ploting the network
plot.mat(M=net, clu=clu)
class(net)<-"mat"
plot(net, clu=clu)
#See corespodning functions for examples for other plotting 
#functions
#presented, that are esentially only the wrappers for "plot.max"
}
\keyword{graphs}% at least one, from doc/KEYWORDS
\keyword{hplot}% at least one, from doc/KEYWORDS
