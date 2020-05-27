#' Sum of Squared deviations from the mean and sum of Absolute Deviations from the median
#' 
#' Functions to compute Sum of Squared deviations from the mean and sum of Absolute Deviations from the median.
#'
#' @param x A numeric vector.
#'
#' @return Sum of Squared deviations from the mean or sum of Absolute Deviations from the median.
#'
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @keywords univar
#' @importFrom stats median
#' 
#' @export

"ss" <-
function(x){sum(x^2)-sum(x)^2/length(x)}

#' @rdname ss
#' 
#' @export

"ad" <-
function(x)sum(abs(x-median(x)))

