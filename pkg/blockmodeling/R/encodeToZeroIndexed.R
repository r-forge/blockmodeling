#' Encode Integer Vector with Consecutive Values Starting from 0
#'
#' This function maps an integer vector to a new vector where the unique values
#' are replaced with consecutive integers starting from 0.
#'
#' @param x An integer vector.
#'
#' @return An integer vector with values from 0 to (number of unique values - 1),
#' maintaining the original order structure.
#'
#' @examples
#' encodeToZeroIndexed(c(10, 20, 10, 30))
#' # Returns: 0 1 0 2
#'
#' @export
encodeToZeroIndexed <- function(x) {
  as.integer(factor(x)) - 1
}
