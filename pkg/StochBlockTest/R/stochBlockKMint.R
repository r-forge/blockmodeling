#' A function for using k-means to initialized the  stochastic one-mode and linked blockmodeling.
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doRNG
#' @import blockmodeling
#' @importFrom stats var
#' @importFrom stats rbinom
#' @importFrom stats kmeans
#'
#' @inheritParams stochBlockORP
#' @param nstart number of random starting points for the classical k-means algorithm (for each set of units). Defaults to \code{100}.
#' @param iter.max Maximum number of iterations from each random starting point (according to \code{nstart}). Defaults to \code{1000}.
#' @param perm Number or partitions obtained by randomly permuting the k-means partition - if 0, no permutations are made, only the original partition is analyzed.
#' @param sharePerm The probability that a unit will have their randomly assigned. Defaults to \code{0.20}.
#'
#' @return A list similar to optRandomParC
#'
#' @export

stochBlockKMint<-function(M,
                          k,
                          nstart=100,
                          iter.max=1000,
                          perm=0,
                          sharePerm=0.20,
                          save.initial.param=TRUE,
                          deleteMs=TRUE,
                          max.iden=10,
                          return.all=FALSE,
                          return.err=TRUE,
                          seed=NULL,
                          RandomSeed=NULL,
                          maxTriesToFindNewPar=perm*10,
                          skip.par=NULL,
                          printRep=ifelse(perm<=10,1,round(perm/10)),
                          n=NULL,
                          nCores=1,
                          useParLapply=FALSE,
                          cl=NULL,
                          stopcl=is.null(cl),
                          ...
){
  dots<-list(...)

  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters


  # if(nCores!=1) {
  #   warning("Mullticore implementation is not yet supported!")
  #   nCores<-1
  # }


  # if(is.null(mingr)){
  #   if(is.null(dots$minUnitsRowCluster)){
  #     mingr<-1
  #   } else {
  #     mingr<-c(dots$minUnitsRowCluster,dots$minUnitsColCluster)
  #   }
  # }
  #
  # if(is.null(maxgr)){
  #   if(is.null(dots$maxUnitsRowCluster)){
  #     maxgr<-Inf
  #   } else {
  #     maxgr<-c(dots$maxUnitsRowCluster,dots$maxUnitsColCluster)
  #   }
  # }

  nmode<-length(k)

  res<-list(NULL)
  err<-NULL
  dots<-list(...)

  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters


  if(is.null(n)) if(nmode==1){
    n<-dim(M)[1]
  } else if(nmode==2){
    n<-dim(M)[1:2]
  } else stop("Number of nodes by modes can not be determined. Parameter 'n' must be supplied!!!")

  if(nmode!=length(n)) stop("The lengths of k and n must match!")

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
                 if(nmode==1){
                   !any(sapply(best.clu,rand2,clu2=res1[[i]]$clu)==1)
                 } else {
                   !any(sapply(best.clu,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(res1[[i]]$clu))==1)
                 }
          )
        ){
          best<-c(best,res1[i])
          best.clu<-c(best.clu,list(res1[[i]]$clu))
        }

        if(length(best)>=max.iden) {
          warning("Only the first ",max.iden," solutions out of ",length(na.omit(err))," solutions with minimal -loglikelihood will be saved.\n")
          break
        }

      }
    }

    names(best)<-paste("best",1:length(best),sep="")

    if(any(na.omit(err)==-Inf) || ss(na.omit(err))!=0 || length(na.omit(err))==1){
      cat("\n\nOptimization of all partitions completed\n")
      cat(length(best),"solution(s) with minimal -loglikelihood =", min(err,na.rm=TRUE), "found.","\n")
    }else {
      cat("\n\nOptimization of all partitions completed\n")
      cat("All",length(na.omit(err)),"solutions have -loglikelihood",err[1],"\n")
    }

    call<-list(call=match.call())
    best<-list(best=best)
    checked.par<-list(checked.par=skip.par)
    if(return.all) res<-list(res=res) else res<-NULL
    if(return.err) err<-list(err=err) else err<-NULL
    if(!exists("initial.param")){
      initial.param<-NULL
    } else initial.param=list(initial.param)

    res<-c(list(M=M),list(ICL=best[[1]][[1]]$ICL),res,best,err,checked.par,call,initial.param=initial.param, list(Random.seed=.Random.seed, cl=cl))
    class(res)<-"opt.more.par"
    return(res)
  })

  dat<-cbind(M, t(M))
  modeVec<-rep(1:nmode, n)
  part<-list()

  for(iMode in 1:nmode){
    iDat<-dat[modeVec==iMode,]
    iDat<-iDat[,colSums(iDat^2)>0]
    iKmRes<-kmeans(iDat, centers = k[iMode],nstart = nstart,iter.max = iter.max)
    part[[iMode]]<-iKmRes$cluster
  }
  perm <- perm+1

  if(nCores==1||!requireNamespace('parallel')){
    if(nCores!=1) {
      oldWarn<-options("warn")
      options(warn=1)
      warning("Only single core is used as package 'parallel' is not available")
      options(warn=oldWarn)
    }


    for(i in 1:perm){
      if(printRep & (i%%printRep==0)) cat("\n\nStarting optimization of the partiton",i,"of",perm,"partitions.\n")
      find.unique.par<-TRUE
      ununiqueParTested=0
      endFun<-FALSE
      while(find.unique.par){
        temppar<-part
        if(i>1){
          for(iPar in 1:length(temppar)){
            while(TRUE){
              tPar<-temppar[[iPar]]
              sel<-rbinom(n[iPar],size=1,prob = sharePerm)
              tPar[sel]<-sample(1:k[iPar],size=sum(sel),replace = TRUE)
              if(length(unique(tPar))==k[iPar]) break
            }
            tPar->temppar[[iPar]]
          }
        } else {
          if(length(temppar)==1)temppar<-temppar[[1]]
          break
        }
        if(length(temppar)==1)temppar<-temppar[[1]]

        find.unique.par<-
          ifelse(is.null(skip.par),
                 FALSE,
                 if(nmode==1) {
                   any(sapply(skip.par,rand2,clu2=temppar)==1)
                 } else any(sapply(skip.par,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(temppar))==1)
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
      res[[i]]<-stochBlock(M=M, clu=temppar,  ...)
      if(deleteMs){
        res[[i]]$M<-NULL
      }
      res[[i]]$best<-NULL

      err[i]<-res[[i]]$err
      if(printRep==1) cat("Final -loglikelihood:",err[i],"\n")
      if(printRep==1) cat("Final partition:   ",unlistPar(res[[i]]$clu),"\n")
    }
  } else {
    oneRep<-function(i,M,n,k,part,...){
      temppar<-part
      if(i>1){
        for(iPar in 1:length(temppar)){
          while(TRUE){
            tPar<-temppar[[iPar]]
            sel<-rbinom(n[iPar],size=1,prob = sharePerm)
            tPar[sel]<-sample(1:k[iPar],size=sum(sel),replace = TRUE)
            if(length(unique(tPar))==k[iPar]) break
          }
          tPar->temppar[[iPar]]
        }
      }
      if(length(temppar)==1)temppar<-temppar[[1]]
      #skip.par<-c(skip.par,list(temppar))

      tres <- try(stochBlock(M=M, clu=temppar,  ...))
      # if(class(tres)=="try-error"){
      if(inherits(tres,what = 'try-error')){
        tres<-list("try-error"=tres, err=Inf, startPart=temppar)
      }
      if(deleteMs){
        tres$M<-NULL
      }
      tres$best<-NULL
      return(list(tres))
    }

    if(!requireNamespace('doParallel')|!requireNamespace('doRNG')) useParLapply<-TRUE

    if(nCores==0){
      nCores<-detectCores()-1
    }
    if(perm==1) perm<-nCores
    pkgName<-utils::packageName()
    if(is.null(pkgName)) pkgName<-utils::packageName(environment(fun.by.blocks))
    if(useParLapply) {
      if(is.null(cl)) cl<-makeCluster(nCores)
      clusterSetRNGStream(cl)
      nC<-nCores
      #clusterExport(cl, varlist = c("kmBlock","kmBlockORP"))
      #clusterExport(cl, varlist = "kmBlock")
      clusterExport(cl, varlist = "pkgName", envir=environment())
      clusterEvalQ(cl, expr={requireNamespace(pkgName,character.only = TRUE)})
      res<-parLapplyLB(cl = cl,1:perm, fun = oneRep, M=M,n=n,k=k,part=part,...)
      if(stopcl) stopCluster(cl)
      res<-lapply(res,function(x)x[[1]])
    } else {
      requireNamespace('doParallel')
      requireNamespace('doRNG')
      if(!getDoParRegistered()|(getDoParWorkers()!=nCores)){
        if(!is.null(cl)) {
          #cl<-makeCluster(nCores)
          registerDoParallel(cl)
        } else registerDoParallel(nCores)
      }
      nC<-getDoParWorkers()

      res<-foreach(i=1:perm,.combine=c, .packages=pkgName) %dorng% oneRep(i=i,M=M,n=n,k=k,part=part,...)
      if(!is.null(cl) & stopcl) {
        registerDoSEQ()
        stopCluster(cl)
      }
    }
    err<-sapply(res,function(x)x$err)
  }
}

