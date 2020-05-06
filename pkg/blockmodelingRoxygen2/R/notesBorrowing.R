#' The notes borrowing network between social-informatics students
#'
#' The data come from a survey conducted in May 1993 on 13 social-informatics students (Hlebec, 1996).
#' The network was constructed from answers to the question, âHow often did you borrow notes from this person?â for each of the fellow students.
#' The respondents indicated the frequency of borrowing by choosing (on a computer) a line of length 1â20, where 1 meant no borrowing. 1 was deducted from all answers, so that 0 now means no borrowing.
#' The data was first used for blockmodeling in Å½iberna (2007).
#' 
#' @format The data set is a valued matrix with 13 rows and columns.
#' 
#' @usage data("notesBorrowing")
#' 
#' @examples 
#' data(notesBorrowing)
#'
#' # Plot the network.
#' # (The function plotMat is from blockmodeling package.)
#' # plotMat(nyt)
#' 
#' @keywords datasets
#' 
#' @references Hlebec, V., (1996). \emph{MetodoloÅ¡ke znaÄilnosti anketnega zbiranja podatkov v analizi omreÅ¾ji: Magistersko delo}. FDV, Ljubljana.
#' 
#' Å½iberna, A. (2007). Generalized blockmodeling of valued networks. \emph{Social Networks}, 29, 105â126. https://doi.org/10.1016/j.socnet.2006.04.002
#'
#' @docType data
#' @name notesBorrowing
NULL