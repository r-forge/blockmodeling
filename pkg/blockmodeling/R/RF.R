#' @encoding UTF-8
#' @title Calculate the value of the Relative Fit function
#'
#' @description The function calculates the value of the Relative Fit function.
#' @param res An object returned by the function \code{optRandomParC}.
#' @param m The number of randomized networks for the estimation of the expected value of a criterion function. It has to be as high as possible. Defaults to 10.
#' @param loops Whether loops are allowed in randomized networks or not, default \code{TRUE}.
#' @return
#' \itemize{
#' \item \code{RF} - The value of the Relative Fit function.
#' \item \code{err} - The value of a criterion function that is used for blockmodeling (for empirical network).
#' \item \code{rand.err} - A vector with the values of the criterion funcion that is used for blockmodeling (for randomized networks).
#' }
#' @details The function randomizes an empirical network to compute the value of the Relative Fit function.
#' The networks are ranomized in such a way that the values on the links are randomly relocated.
#' @examples
#' n <- 8 # If larger, the number of partitions increases 
#' # dramatically as does if we increase the number of clusters
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(3, 5))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#' # Install package blockmodeling and then run the following lines.
#' res <- optRandomParC(M = net, k = 2, rep = 10, approaches = "hom", homFun = "ss", blocks = "com")
#' RF(res = res, m = 100, loops = TRUE)
#' @seealso \code{optRandomParC}
#' @author Marjan Cugmas and Ales Ziberna
#' @references Cugmas, M., Žiberna, A., & Ferligoj, A. (2019). Mechanisms Generating Asymmetric Core-Cohesive Blockmodels. Metodološki Zvezki, 16(1), 17-41.
#' @export
RF <- function(res, m = 10, loops = TRUE){
  errs <- vector(length = m)
  for (i in 1:m){
    randomized <- matrix(sample(res$initial.param$M), nrow = nrow(res$initial.param$M))
    if (loops == FALSE){
      diagonalni <- diag(randomized)[diag(randomized) != 0]
      diag(randomized) <- -1
      randomized[sample(which(randomized == 0), replace = FALSE, size = length(diagonalni))] <- sample(diagonalni)
      diag(randomized) <- 0
    }
    if (err(res) != 0){
      par <- res$initial.param
      names(par) <- gsub(pattern = "dots.", replacement = "", fixed = TRUE, x = names(par))
      par$dots <- NULL
      par$M[,] <- randomized
      utils::capture.output(res.rand <- do.call(optRandomParC, args = par))
      errs[i] <- err(res.rand)
    }
  }
  return(list("RF" = ifelse(err(res) == 0, yes = 1, no = 1 - err(res)/mean(errs)),
              "err" = err(res),
              "rand.err" = unlist(ifelse(err(res) == 0, yes = NA, no = list(errs)))))
}




