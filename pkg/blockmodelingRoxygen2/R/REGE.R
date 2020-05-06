#' REGE - Algorithms for compiting (dis)similarities in terms of regular equivalnece
#' 
#' REGE - Algorithms for compiting (dis)similarities in terms of regular equivalnece (White & Reitz, 1983). 
#' \code{REGE, REGE.for} - Classical REGE or REGGE, as also implemented in Ucinet. Similarities in terms of regular equivalence are computed.  The \code{REGE.for} is a wrapper for calling the FORTRAN subrutine written by White (1985a), modified to be called by R. The \code{REGE} does the same, however it is written in R. The functions with and without ".for" differ only in whether they are implemented  in R of FORTRAN. Needless to say, the functions implemented in FORTRAN are much faster.
#' \code{REGE.ow, REGE.ow.for} - The above function, modified so that a best match is searched for each arc separately (and not for both arcs, if they exist, together). 
#' \code{REGE.nm.for} - REGE or REGGE, modified to use row and column normalized matrices instead of the original matrix.
#' \code{REGE.ownm.for} - The above function, modified so that a best match for an outgoing ties is searched on row-normalized network and for incoming ties on column-normalized network.
#' \code{REGD.for} - REGD or REGDI, a dissimilarity version of the classical REGE or REGGE. Dissimilarities in terms of regular equivalence  are computed.  The \code{REGD.for} is a wrapper for calling the FORTRAN subroutine written by White (1985b), modified to be called by R.
#' \code{REGE.FC}  - Actually an earlier version of REGE. The difference is in the denominator. See Žiberna (2007) for details.
#' \code{REGE.FC.ow} - The above function, modified so that a best match is searched for each arc separately (and not for both arcs, if they exist, together).
#' other - still in testing stage.
#' 
#' @aliases REGE.for REGE.nm.for REGE.ow REGE.ow.for REGE.ownm.for REGD.for REGD.ow.for REGE.FC REGE.FC.ow REGD.ne.for REGD.ow.ne.for REGE.ne.for REGE.nm.diag.for REGE.nm.ne.for REGE.ow.ne.for REGE.ownm.diag.for REGE.ownm.ne.for
#' 
#' 
#' @param M Matrix or a 3 dimensional array representing the network. The third dimension allows for several relations to be analyzed.
#' @param E Initial (dis)similarity in terms of regular equivalnece.
#' @param iter The desired number of iterations.
#' @param until.change Should the iterations be stopped when no change occurs.
#' @param use.diag Should the diagonal be used. If \code{FALSE}, all diagonal elements are set to 0.
#' @param normE Should the equivalence matrix be normalized after each iteration.
#'
#' @return
#'   \item{E}{A matrix of (dis)similarities in terms of regular equivalnece.}
#'   \item{Eall}{An array of (dis)similarity matrices in terms of regular equivalence, each third dimension represets one iteration. For ".for" functions, only the initial and the final (dis)similarities are returned.}
#'   \item{M}{Matrix or a 3 dimensional array representing the network used in the call.}
#'   \item{iter}{The desired number of iterations.}
#'   \item{use.diag}{Should the diagonal be used - for functions implemented in R only.}
#'   ...
#'   
#' @references
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' 
#' White, D. R., & Reitz, K. P. (1983). Graph and semigroup homomorphisms on networks of relations. Social Networks, 5(2), 193-234.
#'
#' White, D. R.(1985a). DOUG WHITE'S REGULAR EQUIVALENCE PROGRAM. Retrieved from \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/REGGE.FOR}
#'
#' White, D. R. (1985b). DOUG WHITE'S REGULAR DISTANCES PROGRAM. Retrieved from \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/REGDI.FOR}
#'
#' White, D. R. (2005). REGGE. Retrieved from \url{http://eclectic.ss.uci.edu/~drwhite/REGGE/}
#' 
#'  #' @author \enc{Aleš Žiberna}{Ales Ziberna} based on Douglas R. White's original REGE and REGD
#' @seealso \code{\link{sedist}}, \code{\link{critFunC}}, \code{\link{optParC}}, \code{\link{plot.mat}}
#'
#' @examples
#' n <- 20
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(5, 15))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- 0
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1) * sample(c(0, 1),
#'    size = tclu[1] * tclu[2], replace = TRUE, prob = c(3/5, 2/5))
#' net[clu == 2, clu == 1] <- 0
#' net[clu == 2, clu == 2] <- 0
#'
#' D <- REGE.for(M = net)$E # Any other REGE function can be used
#' plot.mat(net, clu = cutree(hclust(d = as.dist(1 - D), method = "ward.D"),
#'    k = 2))
#' # REGE returns similarities, which have to be converted to
#' # disimilarities
#'
#' res <- optRandomParC(M = net, k = 2, rep = 10, approaches = "hom", homFun = "ss", blocks = "reg")
#' plot(res) # Hopefully we get the original partition
#' 
#' @keywords cluster graphs 
#' @importFrom stats as.dist
 
"REGE" <-
function(M,E=1,iter=3,until.change=TRUE,use.diag=TRUE){
	n<-dim(M)[1]
	if(n!=dim(M)[2]) stop("M must be a 1-mode matrix")
	if(!use.diag)diag(M)<-0
	Eall<-array(NA,dim=c(n,n,iter+1)) #An array of 'iter' similiaritie matrices
	Eall[,,1]<-E
	diag(Eall[,,1])<-1
	Match<-array(NA,dim=rep(n,4))
	for(i in 2:n){
		for(j in 1:(i-1)){
			for(k in 1:n){
				for(m in 1:n){
					Match[i,j,k,m]<-min(M[i,k],M[j,m]) + min(M[k,i],M[m,j])
					Match[j,i,k,m] <- min(M[j,k],M[i,m]) + min(M[k,j],M[m,i])#/max(1,(max(M[i,k],M[j,m]) + max(M[k,i],M[m,j])+max(M[j,k],M[i,m]) + max(M[k,j],M[m,i])))
				}
			}	
		}
	}

	for(it in 1:iter){
		for(i in 2:n){
			for(j in 1:(i-1)){
				num<-0
				for(k in 1:n){
					#sim<-max(Eall[k,,it]*Match[i,j,k,])
					num<-num+max(Eall[k,,it]*Match[i,j,k,])+max(Eall[k,,it]*Match[j,i,k,])
					#if(i==2&j==1)cat("num = ", num,", den = ",den,", k = ",k,", Maxm1 = ",Maxm1,", ms1 = ",ms1,", Maxm2 = ",Maxm2,", ms2 = ",ms2,"\n")
				}
				#cat("iter=",it,", i=",i,", j=",j,", num=",num,", den=", den,"\n")
				den<-sum(M[c(i,j),])+sum(M[,c(i,j)])
				if(den!=0) {
					Eall[j,i,it+1]<-Eall[i,j,it+1]<-num/den
				} else Eall[j,i,it+1]<-Eall[i,j,it+1]<-1
			}
		}
		diag(Eall[,,it+1])<-1
		if(until.change & all(Eall[,,it]==Eall[,,it+1])){
			Eall<-Eall[,,1:(it+1)]
			break
		}
	}
	itnames<-0:(it)
	itnames[1]<-"initial"
	itnames[it+1]<-"final"
	dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],itnames)
	return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter,use.diag=use.diag))
}

