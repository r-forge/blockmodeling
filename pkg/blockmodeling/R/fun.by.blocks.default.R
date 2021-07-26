#' @rdname funByBlocks
#' @export
"funByBlocks.default" <-
function(x = M, clu, M = x, ignore.diag = "default", sortNames = TRUE, FUN = "mean", ...)
{
    M<-as.array(M)
    dM<-dim(M)
    nn<-ifelse(length(dM)==2,1,dM[1])

    if(is.list(clu)){
      nmode<-length(clu)
      if(nmode>2){
      clu<-unlist(clu)
      clu<-list(clu,clu)
      }
    } else {
      clu<-list(clu,clu)
      nmode<-1
    }
    clu<-lapply(clu,factor)
    if(ignore.diag =="default"){
        if(length(dM)==3){
            ignore.diag <-all(apply(M,1,function(x)identical(ss(diag(x)),0)))&(nmode==1)
        } else ignore.diag <-identical(ss(diag(M)),0)&(nmode==1)  
    }   
    
    if(sortNames) {
        k <- lapply(clu,function(x)sort(unique(x)))
    }else {
        k <- lapply(clu,function(x)unique(x))
    }
    IM.V <- array(NA, dim=c(nn,length(k[[1]]),length(k[[2]])))
    dimnames(IM.V)<-c(list(1:nn),k)
    for(iNet in 1:nn){
        if(length(dM)==3) iM <- M[iNet,,] else iM<-M
        for (i in k[[1]]) {
            for (j in k[[2]]) {
                B<-iM[clu[[1]] == i, clu[[2]] == j, drop = FALSE]
                if (nmode==1 && i == j &&  ignore.diag) diag(B) <- NA
                #removed "dim(B)[1] > 1 &&" from condition above - produces NA's in IM in the diagonal blocks in case of dimension 1x1
                lpar<-list(x = B,...)
                FUNchar<-FUN
                if(!is.character(FUNchar)) FUNchar<-deparse(substitute(FUN))
                if(FUNchar %in% c("mean","sum","min","max")){
                  if(!("na.rm"%in%names(lpar))) lpar<-c(lpar, list(na.rm=TRUE))
                }
                IM.V[iNet,i, j] <- do.call(FUN, lpar)#, na.rm = TRUE
            }
        }
    }
    if(nn==1) return(IM.V[1,,]) else return(IM.V)
}

