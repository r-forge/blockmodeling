#' Computing the threshold
#' 
#' The functions compute the maximum  value of \code{m/cut} where a certain  block is still classified as \code{alt.blocks} and not "null".
#' The difference between \code{find.m} and \code{find.m2} it that \code{find.m} uses an optimization  approach and is faster and more precise
#' than \code{find.m2}. However, \code{find.m} only supports regular ("reg") and complete ("com") as \code{alt.blocks}, while \code{find.m2} supports all block types.
#' Also, \code{find.m} does not always work, especially if \code{cormet} is not "none".
#' 
#' @aliases find.m find.m2 find.cut
#' 
# #' @usage find.m(M, clu, alt.blocks = "reg", diag = !is.list(clu),
# #' cormet = "none", half = TRUE, FUN = "max")
# #' find.m2(M, clu, alt.blocks = "reg", neval = 100, half = TRUE,
# #' ms = NULL, ...)
# #' find.cut(M, clu, alt.blocks = "reg", cuts = "all", ...)
#'
#' @param M A matrix representing the (usually valued) network. For now, only one-relational networks are supported.
#' The network can have one or more modes (different  kinds of units with no ties among themselves.
#' If the network is not two-mode, the matrix must be square.
#' @param clu A partition. Each unique value represents one cluster.
#' If the network  is one-mode, then this should be a vector, else a list of vectors, one for each mode.
#' @param alt.blocks Only one of allowed blocktypes, as alternative to the null block:\cr
#' "com" - complete block\cr
#' "rdo", "cdo" - row and column-dominant blocks (binary, valued, and implicit approach only)\cr
#' "reg" - (f-)regular block\cr
#' "rre", "cre" - row and column-(f-)regular blocks\cr
#' "rfn", "cfn" - row and column-dominant blocks (binary, valued, and implicit approach only)\cr
#' "den" - density block (binary approach only)\cr
#' "avg" - average block (valued approach only).
#' @param diag (default = \code{TRUE}) Should the special status  of diagonal be acknowledged.
#' @param cormet Which method should be used to correct for different maximum  error contributions\cr
#' "none" - no correction\cr
#' "censor" - censor values larger than \code{M}\cr
#' "correct" -  so that the maximum  possible error contribution of the cell is the same regardless of a condition (either that something  must be 0 or at least \code{M}).
#' @param FUN (default = "max") Function f used in row-f-regular, column-f-regular, and f-regular blocks.
#' @param cuts The cuts, which should be evaluated. If \code{cuts="all"} (default), all unique values are evaluated.
#' @param neval A number of different \code{m} values to be evaluated.
#' @param half Should the returned value of m be one half of the value where the inconsistencies are the same.
#' @param ms The values of m where the function should be evaluated.
#' @param \dots Other parameters to \code{crit.fun}.
#'
#' @return A matrix of maximal \code{m/cut} values.
#'
#' @references
#' Doreian, P., Batagelj, V. & Ferligoj, A. \enc{Anuška}{Anuska} (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press.
#'
#' \enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002
#'
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' 
#' @seealso \code{\link{critFunC}} and maybe also \code{\link{optParC}}, \code{\link{plotMat}}
#' 
#' @keywords cluster
#' @importFrom stats optimize

"find.m" <-
function(
	M,	#matrix of a network
	clu,	#partition
	alt.blocks="reg", #alternative block to null block (for now only "reg" is supported)
	diag=!is.list(clu) ,#allow diagonal blocks
	cormet="none", #should we correct for diferent maxismum error contributins
			# "censor" - censor values larger than m
			# "correct" -  so that the maxsimum possible error contribution of the cell is the same regardles of a condition (either that somthing must be o or at least m)
	half = TRUE,	# should the returned value of m be one half of the value where the incosnistencies are the same, otherwise, the m is restricted to max(M)
	FUN="max"
){
	mx<-max(M)*(1+ half)
	mn<-min(M)
  diag=diag
	if(is.list(clu)){
		k<-sapply(clu,function(x)length(unique(x)))
		clu<-lapply(clu,function(x)as.integer(factor(x)))
		if(length(k)>2) {
			for(i in 2:length(clu)){
				clu[[i]]<-clu[[i]] + max(clu[[i-1]])
	  		}
	  		clu<-unlist(clu)
	  		clu<-list(clu,clu)
	  	}
	} else {
		clu<-as.integer(factor(clu))
		clu<-list(clu,clu)
		k<-sapply(clu,function(x)length(unique(x)))
	}
	
	m<-matrix(NA,nrow=k[1],ncol=k[2])
	err<-list(
		reg=function(B,m,FUN){
			nr<-dim(B)[1]	#numer of rows
			nc<-dim(B)[2]	#numer of colums
			sr<-apply(B,1,FUN);er<-m-sr[sr<m]
			sc<-apply(B,2,FUN);ec<-m-sc[sc<m]
			return(sum(er)*nc+ sum(ec)*nr - sum(pmin(rep(er,times=length(ec)),rep(ec,each=length(er)))))	#regular block error
		},
		com=function(B,m,FUN){sumpos(m-B)}
	)
	
	errd<-list(
		com=function(B,m,FUN){sumpos(m-B) + min(0,sum(diag(B))-sumpos(m-diag(B)))},
		reg=function(B,m,FUN){
			nr<-dim(B)[1]	#numer of rows
			nc<-dim(B)[2]	#numer of colums
			sr<-apply(B,1,FUN);er<-m-sr[sr<m]
			sc<-apply(B,2,FUN);ec<-m-sc[sc<m]
			return(sum(er)*nc+ sum(ec)*nr - sum(pmin(rep(er,times=length(ec)),rep(ec,each=length(er))))) #regular block error
		}
	)
	errd.null<-function(B,m){sum(B) - min(0,sum(diag(B))-sumpos(m-diag(B)))}
	
	for(i in 1:k[1]){
		for(j in 1:k[2]){
			B<-M[clu[[1]]==i,clu[[2]]==j, drop=FALSE]
			if(ss(B)==0) m[i,j]<-min(2*B[1,1],mx) else if(i==j&&diag&&sum(dim(B))>1){
				if(errd.null(B,m=mx)>=errd[[alt.blocks]](B,mx,FUN)*ifelse(cormet=="correct",(mx - 0)/(mx - mn),1)){
					m[i,j]<-mx
				}else{ 
					m[i,j]<-optimize(f=function(m,B,alt.blocks,FUN,cormet,mx,mn){corf<-ifelse(cormet=="correct", (mx - 0)/(m - mn),1); if(cormet=="censor") B[B>m]<-m;(errd.null(B,m)-errd[[alt.blocks]](B,m,FUN)*corf)^2},lower=ifelse(cormet=="censor",mn,0),upper=mx,B=B,FUN=FUN,alt.blocks=alt.blocks,cormet=cormet,mx=mx,mn=mn)$minimum
					if(cormet=="correct" && errd.null(B)<err[[alt.blocks]](B,m[i,j],FUN)*(mx - 0)/(m[i,j] - mn)) m[m<=min(B)]<-0				
			}
				
			}else{
				if(sum(B)>=err[[alt.blocks]](B,mx,FUN)*ifelse(cormet=="correct",(mx - 0)/(mx - mn),1)){
					m[i,j]<-mx 
				}else{
					m[i,j]<-optimize(f=function(m,B,alt.blocks,FUN,cormet,mx,mn){corf<-ifelse(cormet=="correct", (mx - 0)/(m - mn),1); if(cormet=="censor") B[B>m]<-m;(sum(B)-err[[alt.blocks]](B,m,FUN)*corf)^2},lower=ifelse(cormet=="censor",mn,0),upper=mx,B=B,FUN=FUN,alt.blocks=alt.blocks,cormet=cormet,mx=mx,mn=mn)$minimum
					if(cormet=="correct" && sum(B)<err[[alt.blocks]](B,m[i,j],FUN)*(mx - 0)/(m[i,j] - mn)) m[m<=min(B)]<-0				
				}
			}
		}
	}
	if(cormet=="censor") m[m<min(M[M>0])]<-0
	if(half) m<-m/2
	return(m)
}

