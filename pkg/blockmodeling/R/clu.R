#' @encoding UTF-8
#' @title Function for extraction of some elements for objects, returend by functions for Generalized blockmodeling
#' 
#' @description
#' Functions for extraction of partition (\code{clu}), all best partitions (\code{partitions}),
#' image or blockmodel (\code{IM})) and  total error or inconsistency (\code{err}) for objects,
#' returned by functions \code{\link{critFunC}} or \code{\link{optRandomParC}}.
#' 
# #' @usage clu(res, which = 1, ...)
# #' partitions(res)
# #' IM(res, which = 1, drop=TRUE, ...)
# #' EM(res, which = 1, drop=TRUE, ...)
# #' err(res, ...) partitions(res)#' 
#' 
#' @param res Result of function \code{\link{critFunC}} or \code{\link{optRandomParC}}.
#' @param which From \code{which} (if there are more than one) "best" solution should the
#' element be extracted. Warning! \code{which} grater than the number of "best" partitions
#' produces an error.
#' @param drop If \code{TRUE} (default), dimensions that have only one level are dropped
#' (\code{drop} function is applied to the final result).
#' @param \dots Not used.
#'
#' 
#' @return The desired element.
#' 
#' 
#' @references
#' Doreian, P., Batagelj, V., & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press.
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' 
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' 
#' @seealso \code{\link{critFunC}}, \code{\link{plot.mat}}, \code{\link{optRandomParC}}
#'
#'
#' @examples
#' n <- 8 # If larger, the number of partitions increases dramatically,
#' # as does if we increase the number of clusters
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(3, 5))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#' 
#' # We select a random partition and then optimize it
#' all.par <- nkpartitions(n = n, k = length(tclu))
#' # Forming the partitions
#' all.par <- lapply(apply(all.par, 1, list),function(x) x[[1]])
#' # to make a list out of the matrix
#' res <- optParC(M = net,
#'    clu = all.par[[sample(1:length(all.par), size = 1)]],
#'     approaches = "hom", homFun = "ss", blocks = "com")
#' plot(res) # Hopefully we get the original partition
#' clu(res) # Hopefully we get the original partition
#' err(res) # Error
#' IM(res) # Image matrix/array.
#' EM(res) # Error matrix/array.
#'
#'  
#' @keywords manip
#' 
#' @export

"clu" <-
function(res,which=1,...){
  if("clu" %in% names(res)){
    res$clu
  }else res$best[[which]]$clu
}

#' @rdname clu
#' 
#' @export

"partitions" <- 
function(res)lapply(res$best,function(x)x$clu)


#' @rdname clu
#' 
#' @export

"err" <-
function(res,...){
  if(is.null(res$best[[1]]$err)){
    min(res$err)
  }else res$best[[1]]$err
}

#' @rdname clu
#' 
#' @export

"IM" <-
function(res,which=1, drop=TRUE, ...){
  if(inherits(res,c("optMorePar","opt.more.par"))){
    IM<-res$best[[which]]$IM
  } else IM<-res$IM
  if(drop)IM<-drop(IM)
  return(IM)
}

#' @rdname clu
#'
#' @export

"EM" <-
function(res,which=1, drop=TRUE,...){
    if(inherits(res,c("optMorePar","opt.more.par"))){
      EM<-res$best[[which]]$EM
    } else EM<-res$EM
    if(drop)EM<-drop(EM)
    return(EM)
}




