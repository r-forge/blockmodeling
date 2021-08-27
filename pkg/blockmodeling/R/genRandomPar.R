#' The function for generating random partitions
#' 
#' The function generates random partitions. The function is meant to be called by the function \code{\link{optRandomParC}.}
#' 
# #' @usage genRandomPar(k, n, seed = NULL, mingr = 1, maxgr = Inf,
# #' addParam = list(genPajekPar = TRUE, probGenMech = NULL))
#'
#' @param k Number of clusters (by modes).
#' @param n Number of units (by modes).
#' @param seed Seed for generating random numbers (partitions).
#' @param mingr Minimal allowed group size.
#' @param maxgr Maximal allowed group size.
#' @param addParam This has to be a list with the following parameters (any or all can be missing, then the default values (see usage) are used):\cr
#' "genPajekPar" - Should the partitions be generated as in Pajek (Batagelj & Mrvar, 2006). If \code{FALSE}, all partitions are selected completely at random while making sure that the partitions have the required number of clusters. \cr
#' \code{probGenMech} - Here the probabilities for 4 different generating mechanisms can be specified. If this is not specified, the value is set to \code{c(1/3, 1/3, 1/3, 0)} if \code{genPajekPar} is \code{TRUE} and to \code{c(0, 0, 0, 1)} if \code{genPajekPar} is \code{FALSE}. The first 3 mechanisms are the same as implemented in Pajek (the second one has almost all units in only one cluster) and the fourth is completely random (from uniform distribution).
#'
#' @return A random partition in the format required by \code{\link{optRandomParC}}. If a network has several modes, then a list of partitions, one for each mode.
#'
#' @references Batagelj, V., & Mrvar, A. (2006). Pajek 1.11. Retrieved from http://vlado.fmf.uni-lj.si/pub/networks/pajek/
#' @author \enc{Aleš Žiberna}{Ales Ziberna}
#' @keywords cluster
#' 
#' @export

"genRandomPar" <-
function(
k,#number of clusters/groups
n,#the number of units in each mode
seed=NULL,#the seed for random generation of partitions
mingr=1,	#minimal alowed group size
maxgr=Inf,	#maximal alowed group size
addParam = list(
	genPajekPar = TRUE, 	#Should the partitions be generated as in Pajek (the other options is completly random)
	probGenMech = NULL) #Here the probabilities for the 4 different mechanizems for specifying the partitions are set. It should be a numeric vector of length 4. If not set this is determined based on the previous parameter.
){
    if(is.null(addParam$probGenMech)){
    	if(is.null(addParam$genPajekPar)||addParam$genPajekPar) probGenMech <- c(1/3,1/3,1/3,0) else probGenMech <- c(0,0,0,1)
    } else probGenMech<-addParam$probGenMech
    if(!is.null(seed))set.seed(seed)
    nmode <- length(k)
    if(nmode==1){
        mingr<-mingr[1]
        maxgr<-maxgr[1]
        find.new.par<-TRUE
		while(find.new.par){
      	  ver<-sample(1:4,size=1,prob=probGenMech)
		  if(k==n) ver<-4
		  if(ver!=4){
			temppar<-integer(n)

			if(ver==1){
				temppar<-1:n%%k+1	
			}

			if(ver==2){
				temppar[1:k]<-1:k
				temppar[(k+1):n]<-k
			}

			if(ver==3){
				temppar[1:k]<-1:k
				temppar[(k+1):n]<-1+trunc(k*runif(n-k))
			}

			for(ii in n:2){
				jj<-trunc(ii*runif(1))
				temppar[c(ii,jj)]<-temppar[c(jj,ii)]
			}
		  }else temppar<-sample(1:k,n,replace=TRUE)
			  temptab<-table(temppar)
			  if((length(temptab)==k)&(min(temptab)>=mingr)&(max(temptab)<=maxgr)){
			  	find.new.par<-FALSE
			  	temppar<-as.numeric(factor(temppar,levels=sample(1:k)))
			  }
		  }
      }else{
        temppar<-NULL
        mingr<-rep(mingr,length.out=nmode)
        maxgr<-rep(maxgr,length.out=nmode)
        for(imode in 1:nmode){
          find.new.par<-TRUE
          while(find.new.par){
          	ver<-sample(1:4,size=1,prob=probGenMech)
            if(ver!=4){
				itemppar<-integer(n[imode])

				if(ver==1){
					itemppar<-1:n[imode]%%k[imode]+1	
				}

				if(ver==2){
					itemppar[1:k[imode]]<-1:k[imode]
					itemppar[(k[imode]+1):n[imode]]<-k[imode]
				}

				if(ver==3){
					itemppar[1:k[imode]]<-1:k[imode]
					itemppar[(k[imode]+1):n[imode]]<-1+trunc(k[imode]*runif(n[imode]-k[imode]))
				}

				for(ii in n[imode]:2){
					jj<-trunc(ii*runif(1))
					itemppar[c(ii,jj)]<-itemppar[c(jj,ii)]
				}
            }else itemppar<-sample(1:k[imode],n[imode],replace=TRUE)
            temptab<-table(itemppar)
            if((length(temptab)==k[imode])&(min(temptab)>=mingr[imode])&(max(temptab)<=maxgr[imode])){
            	find.new.par<-FALSE
            	itemppar<-as.numeric(factor(itemppar,levels=sample(1:k[imode])))
            }
          }
          temppar<-c(temppar,list(itemppar))
        }
      }
      return(temppar)
}

