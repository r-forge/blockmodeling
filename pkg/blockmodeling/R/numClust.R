#' Number of Clusters in a Partition or List of Partitions
#'
#' Computes the number of clusters (unique groups) in a partition or a list of partitions.
#' 
#' If a single vector is provided, it returns the number of unique cluster labels.
#' If a list of partitions is provided (e.g., when dealing with multiple sets of units), 
#' it returns a vector of the number of clusters for each partition in the list.
#'
#' @param clu A vector representing a partition or a list of such vectors (partitions). 
#' Each partition is expected to be a vector of cluster memberships.
#'
#' @return An integer (if `clu` is a single partition) or a vector of integers (if `clu` is a list),
#' representing the number of clusters in each partition.
#'
#' @examples
#' # Single partition
#' numClust(c(1, 1, 2, 2, 3))  # returns 3
#'
#' # List of partitions
#' numClust(list(c(1, 1, 2), c(1, 2, 3, 4)))  # returns c(2, 4)
#'
#' @export
numClust <- function(clu) {
  if (is.list(clu)) {
    return(sapply(clu, numClust))
  } else {
    return(length(unique(clu)))
  }
}