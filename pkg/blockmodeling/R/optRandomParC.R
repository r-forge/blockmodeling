#' @encoding UTF-8
#' @title Optimizing a set of partitions based on the value of a criterion function
#' 
#' @description
#' The function optimizes a set of partitions based on the value of a criterion function  (see \code{\link{critFunC}} for details on the criterion function) for a given network and blockmodel for Generalized blockmodeling (Žiberna, 2007) based on other parameters (see below).
#' The optimization is done through local optimization, where the neighborhood of a partition includes all partitions that can be obtained by moving one unit from one cluster to another or by exchanging two units (from different clusters).
#' A list of paritions can or the number of clusters and a number of partitions to generate can be specified (\code{optParC}).
#'
#' @param k The number of clusters used in the generation of partitions.
#' @param rep The number of repetitions/different starting partitions to check.
#' @param save.initial.param.opt Should the inital parameters(\code{approaches}, ...) of using \code{optParC} be saved. The default value is \code{FALSE}.
#' @param deleteMs Delete networks/matrices from the results of to save space.
#' @param max.iden Maximum number of results that should be saved (in case there are more than \code{max.iden} results with minimal error, only the first \code{max.iden} will be saved).
#' @param switch.names Should partitions that only differ in group names be considered equal. By default it is set to \code{TRUE} if \code{blocks} is either a vector or a list of vectors and to \code{FALSE} otherwise.
#' @param return.all If \code{FALSE}, solution for only the best (one or more) partition/s is/are returned.
#' @param return.err Should the error for each optimized partition be returned.
#' @param seed Optional. The seed for random generation of partitions.
#' @param RandomSeed Optional. Integer vector, containing the random number generator. It is only looked for in the user's workspace.
#' @param parGenFun The function (object) that will generate random partitions. The default function is   \code{\link{genRandomPar}}. The function has to accept the following parameters: \code{k} (number o of partitions by modes, \code{n} (number of units by modes), \code{seed} (seed value for random generation of partition), \code{addParam} (a list of additional parameters).
#' @param mingr Minimal allowed group size.
#' @param maxgr Maximal allowed group size.
#' @param addParam A list of additional parameters for function specified above. In the usage section they are specified for the default function \code{\link{genRandomPar}}.
#' @param maxTriesToFindNewPar The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before - the default value is \code{rep * 1000}.
#' @param skip.par The partitions that are not allowed or were already checked and should therefore be skipped.
#' @param useOptParMultiC For backward compatibility. May be removed soon. See next argument.
#' @param useMulti Which version of local search should be used. Default is currently \code{FALSE}. If \code{FALSE}, first possible all moves in random order and then all possible exchanges in random order are tried. When a move with lower value of criterion function is found, the algorithm moves to this new partition. If \code{TRUE} the version of local search where all possible moves and exchanges are tried first and then the one with the lowest error is selected and used. In this case, several optimal partitions are found. \code{maxPar} best partitions are returned.
#' @param printRep Should some information about each optimization be printed.
#' @param n The number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is one mode (both modes are in rows and columns).
#' @param nCores Number of cores to be used. Value \code{0} means all available cores. It can also be a cluster object.
#' @param useParLapply Should \code{parLapplyLB} be used (otherwise \code{mforeach} is used). Defaults to true as it needs less dependencies. It might be removed in future releses and only allow the use of parLapplyLB.
#' @param  cl The cluster to use (if formed beforehand). Defaults to \code{NULL}.
#' @param  stopcl Should the cluster be stoped after the function finishes. Defaults to \code{is.null(cl)}.
#' @param genPajekPar Should the partitions be generated as in Pajek.
#' @param probGenMech Should the probabilities for different mechanisms for specifying the partitions be set. If \code{probGenMech} is not set, it is determined based on the parameter \code{genPajekPar}.
#' @param \dots Arguments passed to other functions, see \code{\link{critFunC}}.
#' @inheritParams critFunC 
#'
#' @return
#'   \item{M}{The matrix of the network analyzed.}
#'   \item{res}{If \code{return.all = TRUE} - A list of results the same as \code{best} - one \code{best} for each partition optimized.}
#'   \item{best}{A list of results from \code{optParC}, only without \code{M}.}
#'   \item{err}{If \code{return.err = TRUE} - The vector of errors or inconsistencies of the empirical  network with the ideal network for a given blockmodel (model,approach,...) and parititions.}
#'   \item{nIter}{The vector of the number of iterations used - one value for each starting partition that was optimized. It can show that \code{maxiter} is too low if a lot of these values have the value of \code{maxiter}.}
#'   \item{checked.par}{If selected - A list of checked partitions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includes the partitions in \code{skip.par}.}
#'   \item{call}{The call used to call the function.}
#'   \item{initial.param}{If selected - The initial parameters are used.}
#'   
#' @section Warning:
#' It should be noted that the time complexity of package blockmodeling is increasing with
#' the number of units and the number of clusters (due to its algorithm). Therefore the analysis
#' of network with more than 100 units can take a lot of time (from a few hours to a few days).
#' 
#' @references Batagelj, V., & Mrvar, A. (2006). Pajek 1.11. Retrieved from \url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}
#' 
#' Doreian, P., Batagelj, V. & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press.
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002
#' 
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' @seealso \code{\link{critFunC}}
#' 
#' @examples
#' n <- 8 # If larger, the number of partitions increases dramatically
#' # as does if we increase the number of clusters
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(3, 5))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#'
#' # Optimizing 10 random chosen partitions with optRandomParC
#' res <- optRandomParC(M = net, k = 2, rep = 10,
#' approaches = "hom", homFun = "ss", blocks = "com")
#' plot(res) # Hopefully we get the original partition
#' 
#' @keywords cluster graphs
#' @import methods
#' @import parallel
#' @importFrom stats na.omit runif
#'
#' @export

"optRandomParC" <-function(M, 
                           k,#number of clusters/groups
                           approaches, #generalized blockmodeling approach
                           blocks, #allowed block types as a vector, list or array.
                           rep,#number of repetitions/different starting partitions to check
                           save.initial.param=TRUE,  #save the initial parametrs of this call
                           save.initial.param.opt=FALSE,  #save the initial parametrs for calls to optParC
                           deleteMs=TRUE, #delete networks/matrices from results of optParC
                           max.iden=10, #the maximum number of results that should be saved (in case there are more than max.iden results with minimal error, only the first max.iden will be saved)
                           switch.names=NULL,#should partitions that only differ in group names be considert equal (is c(1,1,2)==c(2,2,1))
                           return.all=FALSE,#if 'FALSE', solution for only the best (one or more) partition/s is/are returned
                           return.err=TRUE,#if 'FALSE', only the resoults of critFun are returned (a list of all (best) soulutions including errors), else the resoult is list
                           seed=NULL,#the seed for random generation of partitions
                           RandomSeed=NULL, # the state of .Random.seed (e.g. as saved previously). Should not be "typed" by the user
                           parGenFun = genRandomPar, #The function that will generate random partitions. It should accept argumetns: k (number of partitions by modes, n (number of units by modes), seed (seed value for random generation of partition), addParam (a list of additional parametres)
                           mingr=NULL, #minimal alowed group size (defaults to c(minUnitsRowCluster,minUnitsColCluster) if set, else to 1) - only used for parGenFun function 
                           maxgr=NULL, #maximal alowed group size (default to c(maxUnitsRowCluster,maxUnitsColCluster) if set, else to Inf) - only used for parGenFun function 
                           addParam=list(  #list of additional parameters for gerenrating partitions. Here they are specified for dthe default function "genRandomPar"
                             genPajekPar = TRUE,     #Should the partitions be generated as in Pajek (the other options is completly random)
                             probGenMech = NULL),    #Here the probabilities for different mechanizems for specifying the partitions are set. If not set this is determined based on the previous parameter.
                           maxTriesToFindNewPar=rep*10,    #The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before 
                           skip.par = NULL, #partitions to be skiped
                           useOptParMultiC = FALSE, # For backward compatibility. May be removed soon. See next argumetent.
                           useMulti = useOptParMultiC, #Should the "Multi" vesrsion of the optParC functions be used? Defaults to FALSE, which is usually faster, but in a sense not so thorough.
                           printRep= ifelse(rep<=10,1,round(rep/10)), #should some information about each optimization be printed
                           n=NULL, #the number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is onemode (both modes are in rows and columns)
                           nCores=1, #number of cores to be used 0 -means all available cores, can also be a cluster object
						   useParLapply = TRUE, # Should parLapplyLB be used (otherwise foreach is used)
							cl = NULL, #the cluster to use (if formed beforehand)
							stopcl = is.null(cl), # should the cluster be stoped						   
                           ... #paramters to optParC
){
  dots<-list(...) #this might not be need - can be removed and all latter occurencies given sufficent testing. Left for now as there is not enought time.
  if(is.null(switch.names)){
    if(is.list(blocks)){
      switch.names<-all(sapply(blocks,is.vector))
    } else switch.names<-is.vector(blocks)
  }
  
  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
  
  if(is.null(mingr)){
    if(is.null(dots$minUnitsRowCluster)){
      mingr<-1
    } else {
      mingr<-c(dots$minUnitsRowCluster,dots$minUnitsColCluster)
    }
  }
  
  if(is.null(maxgr)){
    if(is.null(dots$maxUnitsRowCluster)){
      maxgr<-Inf
    } else {
      maxgr<-c(dots$maxUnitsRowCluster,dots$maxUnitsColCluster)
    }
  }
  
  nmode<-length(k)
  
  res<-list(NULL)
  err<-NULL
  nIter<-NULL
  
  if(is.null(n)) if(nmode==1){
    n<-dim(M)[1]
  } else if(nmode==2){
    n<-dim(M)[1:2]
  } else warning("Number of nodes by modes can not be determined. Parameter 'n' must be supplied!!!")
  
  if(!is.null(RandomSeed)){
    .Random.seed <-  RandomSeed
  } else if(!is.null(seed))set.seed(seed)
  
  on.exit({
    res1 <- res[which(err==min(err, na.rm = TRUE))]
    best<-NULL
    best.clu<-NULL
    for(i in 1:length(res1)){
      for(j in 1:length(res1[[i]]$best)){
        if(
          ifelse(is.null(best.clu),
                 TRUE,
                 if(nmode==1) ifelse(switch.names,
                                     !any(sapply(best.clu,rand2,clu2=res1[[i]]$clu)==1),
                                     !any(sapply(best.clu,function(x)all(x==res1[[i]]$clu)))
                 ) else ifelse(switch.names,
                               !any(sapply(best.clu,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(res1[[i]]$clu))==1),
                               !any(sapply(best.clu,function(x)all(unlist(x)==unlist(res1[[i]]$clu))))
                 )
          )
        ){
          best<-c(best,res1[i])
          best.clu<-c(best.clu,list(res1[[i]]$clu))
        }
        
        if(length(best)>=max.iden) {
          warning("Only the first ",max.iden," solutions out of ",length(na.omit(err))," solutions with minimal error will be saved.\n")
          break
        }
        
      }
    }
    
    names(best)<-paste("best",1:length(best),sep="")
    
    if(any(na.omit(err)==Inf) || ss(na.omit(err))!=0 || length(na.omit(err))==1){
      cat("\n\nOptimization of all partitions completed\n")
      cat(length(best),"solution(s) with minimal error =", min(err,na.rm=TRUE), "found.","\n")
    }else {
      cat("\n\nOptimization of all partitions completed\n")
      cat("All",length(na.omit(err)),"solutions have err",err[1],"\n")
    }
    
    call<-list(call=match.call())
    best<-list(best=best)
    checked.par<-list(checked.par=skip.par)
    if(return.all) res<-list(res=res) else res<-NULL
    if(return.err) err<-list(err=err) else err<-NULL
    if(!exists("initial.param")){
      initial.param<-NULL
    } else initial.param=list(initial.param)
    
    res<-c(list(M=M),res,best,err,list(nIter=nIter),checked.par,call,initial.param=initial.param, list(Random.seed=.Random.seed))
    class(res)<-"optMorePar"
    return(res)
  })
  
   if(nCores==1||!requireNamespace("parallel")){
     if(nCores!=1) {
       oldWarn<-options("warn")
       options(warn=1)
       warning("Only single core is used as package 'parallel' is not available")
       options(warn=oldWarn)
     }
    for(i in 1:rep){
      if(printRep & (i%%printRep==0)) cat("\n\nStarting optimization of the partiton",i,"of",rep,"partitions.\n")
      find.unique.par<-TRUE
      ununiqueParTested=0
      while(find.unique.par){
        temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
        
        find.unique.par<-
          ifelse(is.null(skip.par),
                 FALSE,
                 if(nmode==1) ifelse(switch.names,
                                     any(sapply(skip.par,rand2,clu2=temppar)==1),
                                     any(sapply(skip.par,function(x)all(x==temppar)))
                 ) else ifelse(switch.names,
                               any(sapply(skip.par,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(temppar))==1),
                               any(sapply(skip.par,function(x)all(unlist(x)==unlist(temppar))))
                 )
          )
        ununiqueParTested=ununiqueParTested+1
        endFun<-ununiqueParTested>=maxTriesToFindNewPar
        if(endFun) {
          break
        } else if(ununiqueParTested%%10==0) cat(ununiqueParTested,"partitions tested for unique partition\n")
      }
      
      if(endFun) break
      
      skip.par<-c(skip.par,list(temppar))
      
      if(printRep==1) cat("Starting partition:",unlistPar(temppar),"\n")
      #if(useOptParMultiC){
      #    res[[i]]<-optParMultiC(M=M, clu=temppar,  save.initial.param= save.initial.param.opt, ...)
      #}else  res[[i]]<-optParC(M=M, clu=temppar,  save.initial.param= save.initial.param.opt,  ...)
      res[[i]]<-optParC(M=M, clu=temppar, approaches=approaches, blocks=blocks, useMulti=useMulti, save.initial.param= save.initial.param.opt,  ...)
      if(deleteMs){
        res[[i]]$M<-NULL
        res[[i]]$resC$M<-NULL
      }
      
      err[i]<-res[[i]]$err
      nIter[i]<-res[[i]]$resC$nIter
      if(printRep==1) cat("Final error:",err[i],"\n")
      if(printRep==1) cat("Final partition:   ",unlistPar(res[[i]]$clu),"\n")
    }
  } else {
    oneRep<-function(i,M,approaches, blocks, n,k,mingr,maxgr,addParam,rep,...){
      if(printRep) cat("\n\nStarting optimization of the partiton",i,"of",rep,"partitions.\n")
      temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
      
      #skip.par<-c(skip.par,list(temppar))
      
      #if(useOptParMultiC){
      #    tres <- try(optParMultiC(M=M, clu=temppar,  save.initial.param= save.initial.param.opt,  ...))
      #}else  tres <- try(optParC(M=M, clu=temppar,  save.initial.param= save.initial.param.opt,  ...))
      tres <- try(optParC(M=M, clu=temppar, approaches=approaches, blocks=blocks, useMulti=useMulti, save.initial.param= save.initial.param.opt,  ...))
      
      if(class(tres)=="try-error"){
        tres<-list("try-error"=tres, err=Inf, nIter=Inf, startPart=temppar)
      }
      if(deleteMs){
        tres$M<-NULL
        tres$resC$M<-NULL
      }
      #            err[i]<-res[[i]]$err
      #            nIter[i]<-res[[i]]$resC$nIter
      return(list(tres))
    }  
     
	if(useParLapply||!requireNamespace("doParallel")||!requireNamespace("foreach")||!requireNamespace("doRNG")) useParLapply<-TRUE 
    if(nCores==0){
       nCores<-detectCores()-1                    
    }
	
	pkgName<-utils::packageName()
	if(is.null(pkgName)) {
		pkgName<-utils::packageName(environment(optParC))
		cat("Package name set by a trick!\n")
	}
	
	if(useParLapply) {
       if(is.null(cl)) cl<-makeCluster(nCores)
       clusterSetRNGStream(cl)
       nC<-nCores
       #clusterExport(cl, varlist = c("kmBlock","kmBlockORP"))
       #clusterExport(cl, varlist = "kmBlock")
       clusterExport(cl, varlist = "pkgName", envir=environment())	   
       clusterEvalQ(cl, expr={require(pkgName,character.only = TRUE)})
       res<-parLapplyLB(cl = cl,1:rep, fun = oneRep,M=M,approaches=approaches, blocks=blocks ,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep,...)
       if(stopcl) stopCluster(cl)
       res<-lapply(res,function(x)x[[1]])
       err<-sapply(res,function(x)x$err)
	   nIter<-sapply(res,function(x)x$resC$nIter)
    } else {
		requireNamespace("doParallel")
		requireNamespace("doRNG")
		requireNamespace("foreach")
		`%dorng%`<-doRNG::`%dorng%`
		`%dopar%`<-foreach::`%dopar%`
		if(!foreach::getDoParRegistered()){
			doParallel::registerDoParallel(nCores)
		}
		nC<-foreach::getDoParWorkers()

		pkgName<-utils::packageName()
		res<-foreach::foreach(i=1:rep,.combine=c, .packages=pkgName) %dorng% oneRep(i=i,M=M,approaches=approaches, blocks=blocks ,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep,...)
		err<-sapply(res,function(x)x$err)
		nIter<-sapply(res,function(x)x$resC$nIter)
	}
  }
}


unlistPar<-function(part){
  if(is.list(part)){
    part<-sapply(part,paste,collapse=" ")
    part<-paste(paste("\nMode ", 1:length(part),":",sep=""), part,collapse="",sep=" ")
  }
  part
}


parArr2VecC<-function(parArr,nUnitsClu=NULL){
  if(is.null(nUnitsClu)){
    nUnitsClu<-apply(parArr,2,function(x)sum(!is.na(x)))
  }
  n<-sum(nUnitsClu)
  nClus <- length(nUnitsClu)
  if(!is.integer(parArr)){
    parArr<-apply(parArr,2,as.integer)
  }
  resC<-.C("parArr2Vec",n=as.integer(n), nClus = nClus, nUnitsClu=as.integer(nUnitsClu),parArr=parArr, parVec=integer(n), NAOK=TRUE)
  return(resC$parVec)
}


parVec2ArrC<-function(parVec){
  n<-length(parVec)
  parVec<-as.integer(as.factor(parVec))- as.integer(1)
  nClus <- as.integer(max(parVec)+1)
  nUnitsClu<-integer(nClus)
  parArr<-matrix(NA,ncol=nClus,nrow=n)
  parArr<-apply(parArr,2,as.integer)
  resC<-.C("parVec2Arr",n=n, nClus = nClus, nUnitsClu=nUnitsClu,parArr=parArr, parVec=parVec, NAOK=TRUE)
  parArr<-resC$parArr
  return(parArr)
}




