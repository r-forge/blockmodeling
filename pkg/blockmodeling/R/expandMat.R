#' @encoding UTF-8
#' @title Expands a square matrix by repeating each row/column the specified number of times.
#' @param mat A square matrix to be exapanded
#' @param nn A vector of number of times each row/column must be repeated. Its length must match the number of rows/columns
#' @return Sum of squared deviations from the mean using only valid (non NA) values.
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' @keywords manip
#' @export
expandMat<-function(mat, nn){
v<-rep(1:length(nn),nn)
mat[v,v]
}
