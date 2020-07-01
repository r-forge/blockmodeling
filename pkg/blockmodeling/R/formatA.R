#' @encoding UTF-8
#' @title A formating function for numbers
#' 
#' @description
#' Formats a vector or matrix of numbers so that all have equal length (digits). This is especially suitable for printing tables.
#'
# #' @usage formatA(x, digits = 2, FUN = round, ...)
#' 
#' @param x A numerical vector or matrix.
#' @param digits The number of desired digits.
#' @param FUN Function used for "shortening" the numbers.
#' @param \dots Additional arguments to \code{format}.
#'
#' @return A character vector or matrix.
#'
#' @examples
#' A <- matrix(c(1, 1.02002, 0.2, 10.3), ncol = 2)
#' formatA(A)
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @seealso \code{\link{find.m}}, \code{\link{find.m2}}, \code{\link{find.cut}}
#' @keywords character
#' 
#' @export

"formatA" <-
function(x,digits=2, FUN=round,...){
	noquote(format(FUN(x, digits=digits),...))
}

