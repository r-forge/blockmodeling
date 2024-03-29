#' @encoding UTF-8
#' @title Computes distances in terms of Structural equivalence (Lorrain & White, 1971)
#' 
#' @description
#' The functions compute the distances in terms of Structural equivalence (Lorrain and White, 1971) between the units of a one-mode network.
#' Several options for treating the diagonal values are supported.
#' 
# #' @usage
# #' sedist(M, method = "default", fun = "default",
# #' fun.on.rows = "default", handle.interaction = "switch",
# #' use = "pairwise.complete.obs", ...)
#'
#' @param M A matrix representing the (usually valued) network. For now, only one-relational networks are supported. The network must be one-mode.
#' @param method The method used to compute distances - any of the methods allowed by functions dist, \code{"cor"} or \code{"cov"} (all \code{package::stats}) or just \code{"cor"} or \code{"cov"} (given as a character).
#' @param fun Which function should be used to compute distances (given as a character).
#' @param fun.on.rows For non-standard function - does the function compute measure on rows (such as \code{"cor"}, \code{"cov"},...) of the data matrix (as opposed to computing measure on columns (such as \code{dist}).
#' @param handle.interaction How should the interaction between the vertices analysed be handled:\cr
#'        \code{"switch"} (the default) - assumes that when comparing units i and j, M[i,i] should be compared with M[j,j] and M[i,j] with M[j,i]. These two comparisons are weighted by 2. This should be used with Euclidean distance to get the corrected Euclidean distance with p = 2.\cr
#'        \code{"switch2"} - the same (alias)\cr
#'        \code{"switch1"} - the same as above, only that the two comparisons are weighted by 1. This should be used with Euclidean distance to get the corrected Wuclidean distance with p = 1.\cr
#'        \code{"ignore"} (diagonal) - Diagonal is ignored. This should be used with Euclidean distance to get the corrected Euclidean distance with p = 0.\cr
#'        \code{"none"} - the matrix is used "as is"
#' @param use For use with methods \code{"cor"} and \code{"cov"}, for other methods (the default option should be used if \code{handle.interaction == "ignore"}), \code{"pairwise.complete.obs"} are always used, if \code{stats.dist.cor.cov = TRUE}.
#' @param \dots Additional arguments to \code{fun}
#' 
#' @details
#' If both \code{method} and \code{fun} are \code{"default"}, the Euclidean distances are computed.
#' The \code{"default"} method for \code{fun = "dist"} is "euclidean" and for \code{fun  = "cor"} "pearson".  
#'
#' @return A matrix (usually of class dist) is returned.
#' 
#' @references
#' Batagelj, V., Ferligoj, A., & Doreian, P. (1992). Direct and indirect methods for structural equivalence. Social Networks, 14(1-2), 63-90. doi: 10.1016/0378-8733(92)90014-X
#' 
#' Lorrain, F., & White, H. C. (1971). Structural equivalence of individuals in social networks. Journal of Mathematical Sociology, 1(1), 49-80. doi: 10.1080/0022250X.1971.9989788
#' 
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @seealso \code{\link{dist}}, \code{\link{hclust}}, \code{\link{REGE}}, \code{\link{optParC}}, \code{\link{optParC}}, \code{\link{optRandomParC}}
#' 
#' @examples
#' # Generating a simple network corresponding to the simple Sum of squares
#' # Structural equivalence with blockmodel:
#' # null com
#' # null null
#' n <- 20
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(5, 15))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#'
#' D <- sedist(M = net)
#' plot.mat(net, clu = cutree(hclust(d = D, method = "ward"), k = 2))
#' 
#' @keywords cluster graphs
#' @importFrom stats as.dist cor cov na.omit
#' 
#' @export

"sedist" <-
function(
	M,	#matrix (of a network)
	method="default", 	# the a method used to compute distances - any of the methods alloed by functions dist, cor or cov {all package::stats} or just "cor" or "cov" (given as character)
	fun="default",	#which function should be used to comput distacnes (given as character),
	fun.on.rows="default", # for non-standard function - does it compute measure on rows (such as cor, cov,...) of the data matrix.
#	stats.dist.cor.cov=TRUE,	#call "stats::dist", "stats::cor" or "stats::cov", not "dist", "cor" or "cov", if nonstandard functions are used, they should exemp the same arguments as those in package stats
	handle.interaction="switch",	#how should the interaction between the vertices analysed be handled:
						# "switch" (the default) - assumes that when comparing units i and j, M[i,i] should be compared with M[j,j] and M[i,j] with M[j,i]
						# "switch1" - the same as above, only that each pair occurs only once
						# "switch2" - an alias for switch
						# "ignore" (diagonal) - Diagonal is ignored
						# "none" - the matrix is used "as is"
	use = "pairwise.complete.obs",	#for use with methods "cor" and "cov", for other methods (the default option should be used if handle.interaction=="ignore"), "pairwise.complete.obs" are always used, if stats.dist.cor.cov=TRUE
	#p=2	,#The power of the Minkowski distance in functin dist if stats.dist.cor.cov=TRUE
	... #other argumets passed to fun
)
{

	method<-match.arg(method, choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski","pearson", "kendall", "spearman","dist","cor", "cov", "default"))
	handle.interaction<-match.arg(handle.interaction, choices=c("switch", "switch1", "switch2", "ignore", "none"))
	if(handle.interaction=="switch2")handle.interaction<-"switch"
	if(any(method=="default", fun=="default")){
		if(all(method=="default", fun=="default")){
			fun<-"dist"
			method<-"euclidean"
		} else if(fun=="default"){
			if(method %in% c("pearson", "kendall", "spearman")) fun<-"cor"
			if(method %in% c("cor", "cov")){
				fun<-method
				method<-"pearson"
			}
			if(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) fun<-"dist"
		} else {
			if(fun %in% c("cor","cov")) method<-"pearson"
			if(fun=="dist") method<-"euclidean"
		}
	}

	if(handle.interaction=="ignore"&& fun %in% c("cor","cov") && use != "pairwise.complete.obs")warning("The option use='pairwise.complete.obs' should be used with handle.interaction=='ignore' && fun %in% c('cor','cov')")

#	if(fun %in% c("dist", "cor" or "cov") && stats.dist.cor.cov) fun<-paste("stats::",fun,sep="")
	if(fun.on.rows=="default") if(fun %in% c("cor","cov")){
		fun.on.rows<-TRUE
	} else fun.on.rows<-FALSE

	n<-dim(M)[1]
	if(n!=dim(M)[2]) stop("This function is suited for one-mode networks only")
    if(fun %in% c("cor", "cov")) usearg<-list(use=use) else usearg<-NULL #usearg

	if(handle.interaction %in% c("switch","switch1")){
		if(fun=="cor"){
			cor1<-function(...)cor(...)[1,2]
			fun<-"cor1"
		}
		if(fun=="cov"){
			cor1<-function(...)cov(...)[1,2]
			fun<-"cov1"
		}
		X<-cbind(M,t(M))
		res<-matrix(NA,ncol=n,nrow=n)
		for(i in 2:n)for(j in seq(length=(i-1))){
			jind<-seq(length=2*n)
			jind[i]<-j
			jind[j]<-i
			jind[n+i]<-ifelse(handle.interaction=="switch",n+j,NA)
			jind[n+j]<-ifelse(handle.interaction=="switch",n+i,NA)
			Xij<-rbind(X[i,],X[j,jind])
			if(fun.on.rows)Xij<-t(Xij)
			res[i,j]<-do.call(fun,args=c(list(x=Xij, method=method,...),usearg))
		}
		if(handle.interaction=="switch1" & fun=="dist" & !(method%in%c("maximum","binary"))) res<-res*sqrt((n-1)/n)
		res<-as.dist(res)
	}else{
		if(handle.interaction=="ignore") diag(M)<-NA
		X<-cbind(M,t(M))
		if(fun.on.rows)X<-t(X)
		res<-do.call(fun,args=c(list(x=X, method=method,...),usearg))
	}
	if(inherits(res,"dist"))attr(res,"Labels")<-rownames(M)
	if(is.matrix(res))dimnames(res)<-dimnames(M)
	return(res)	
}



"sedistX" <-    function(
    X,	#a matrix composed of network and network transposed
    method="default", 	# the a method used to compute distances - any of the methods alloed by functions dist, cor or cov {all package::stats} or just "cor" or "cov" (given as character)
    fun="default",	#which function should be used to comput distacnes (given as character),
    fun.on.rows="default", # for non-standard function - does it compute measure on rows (such as cor, cov,...) of the data matrix.
    #	stats.dist.cor.cov=TRUE,	#call "stats::dist", "stats::cor" or "stats::cov", not "dist", "cor" or "cov", if nonstandard functions are used, they should exemp the same arguments as those in package stats
    handle.interaction="switch",	#how should the interaction between the vertices analysed be handled:
    # "switch" (the default) - assumes that when comparing units i and j, M[i,i] should be compared with M[j,j] and M[i,j] with M[j,i]
    # "switch1" - the same as above, only that each pair occurs only once
    # "switch2" - an alias for switch
    # "ignore" (diagonal) - Diagonal is ignored
    # "none" - the matrix is used "as is"
    use = "pairwise.complete.obs",	#for use with methods "cor" and "cov", for other methods (the default option should be used if handle.interaction=="ignore"), "pairwise.complete.obs" are always used, if stats.dist.cor.cov=TRUE
    #p=2	,#The power of the Minkowski distance in functin dist if stats.dist.cor.cov=TRUE
    ... #other argumets passed to fun
){
    
    method<-match.arg(method, choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski","pearson", "kendall", "spearman","dist","cor", "cov", "default"))
    handle.interaction<-match.arg(handle.interaction, choices=c("switch", "switch1", "switch2", "ignore", "none"))
    if(handle.interaction=="switch2")handle.interaction<-"switch"
    if(any(method=="default", fun=="default")){
        if(all(method=="default", fun=="default")){
            fun<-"dist"
            method<-"euclidean"
        } else if(fun=="default"){
            if(method %in% c("pearson", "kendall", "spearman")) fun<-"cor"
            if(method %in% c("cor", "cov")){
                fun<-method
                method<-"pearson"
            }
            if(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) fun<-"dist"
        } else {
            if(fun %in% c("cor","cov")) method<-"pearson"
            if(fun=="dist") method<-"euclidean"
        }
    }
    
    if(handle.interaction=="ignore"&& fun %in% c("cor","cov") && use != "pairwise.complete.obs")warning("The option use='pairwise.complete.obs' should be used with handle.interaction=='ignore' && fun %in% c('cor','cov')")
    
    #	if(fun %in% c("dist", "cor" or "cov") && stats.dist.cor.cov) fun<-paste("stats::",fun,sep="")
    if(fun.on.rows=="default") if(fun %in% c("cor","cov")){
        fun.on.rows<-TRUE
    } else fun.on.rows<-FALSE
    
    n<-dim(X)[1]
    if(dim(X)[2]%%n!=0) stop("The columns must be a multiple of the rows")
    k<-dim(X)[2]/n
    if(fun %in% c("cor", "cov")) usearg<-list(use=use) else usearg<-NULL #usearg
    
    if(handle.interaction %in% c("switch","switch1")){
        if(fun=="cor"){
            cor1<-function(...)cor(...)[1,2]
            fun<-"cor1"
        }
        if(fun=="cov"){
            cor1<-function(...)cov(...)[1,2]
            fun<-"cov1"
        }
        res<-matrix(NA,ncol=n,nrow=n)
        for(i in 2:n)for(j in seq(length=(i-1))){
            jind<-seq(length=k*n)
            for(l in seq(0,k-1,by = 2)){
                jind[l*n+i]<-j
                jind[l*n+j]<-i
                if((l+1)<k){
                    if(handle.interaction=="switch"){
                        jind[(l+1)*n+i]<-(l+1)*n+j
                        jind[(l+1)*n+j]<-(l+1)*n+i
                    }else{
                        jind[(l+1)*n+i]<-NA
                        jind[(l+1)*n+j]<-NA
                    }
                } 
            }
            Xij<-rbind(X[i,],X[j,jind])
            if(fun.on.rows)Xij<-t(Xij)
            res[i,j]<-do.call(fun,args=c(list(x=Xij, method=method,...),usearg))
        }
        if(handle.interaction=="switch1" & fun=="dist" & !(method%in%c("maximum","binary"))) res<-res*sqrt((n-1)/n)
        res<-as.dist(res)
    }else{
        for(i in 1:n){
            for(l in 0:(k-1)){
                X[i,l*n+i]<-NA
            }
        }
        if(fun.on.rows)X<-t(X)
        res<-do.call(fun,args=c(list(x=X, method=method,...),usearg))
        if(fun.on.rows)X<-t(X)
    }
    if(inherits(res,"dist"))attr(res,"Labels")<-rownames(X)
    if(is.matrix(res))colnames(res)<-rownames(res)<-rownames(X)
    return(res)	
}


