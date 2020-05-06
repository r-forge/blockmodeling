#' Comparing partitions
#' 
#' Rand Index and Rand Index corrected/adjusted for chance for comparing partitions (Hubert & Arabie, 1985).
#' The names of the clusters do not matter.
#' 
#' @aliases crand crand2 rand rand2
#' 
# #' @usage
# #' rand(tab)
# #' rand2(clu1, clu2)
# #' crand(tab)
# #' crand2(clu1, clu2)
#' 
#' @param clu1 The two partitions to be compared, given in the form of vectors, where for each unit a cluster membership is given.
#' @param clu2 The two partitions to be compared, given in the form of vectors, where for each unit a cluster membership is given.
#' @param tab A contingency table obtained as a table(clu1, clu2).
#'
#' @return The value of Rand Index (corrected/adjusted for chance)
#'
#' @references Hubert, L., & Arabie, P. (1985). Comparing Partitions. Journal of Classification, 2(1), 193-218.
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @keywords cluster
#' 

"rand" <-
function (tab) #Hubert & Arabie
{
    n <- sum(tab)
	1 + (sum(tab^2) - (sum(rowSums(tab)^2) + sum(colSums(tab)^2))/2)/choose(n, 2)
}

