#' A function that save an object given to x in a file named with the name of the object and an edning ".Rdata". A different filename can be aslo suplied, but then this is the same as save.
#'
#' @param x An object to be saved.
#' @param fileName A filename into which the object should be saved. Usually omitted and constructed automatically.
saveObj<-function(x,fileName=NULL){
    objName <- deparse(substitute(x))
    if (is.null(fileName) & (!is.null(objName))) 
    objName<-deparse(substitute(x)) 
    if(is.null(fileName)&(!is.null(objName))) fileName<-paste(objName,".Rdata",sep="")
    if(!is.null(fileName)) save(list=objName, file = fileName) else stop("Suitable filename was not provided and could not be created based on object name")
}


#' A function that loads an object with the name used as the first arguemnt from a file named with the name of the object and an edning ".Rdata". A different filename can be aslo suplied, but then this is the same as load. Reverse of saveObj.
#'
#' @param x An object to be loaded (acctually unquoted filename without the ".Rdata" ending from which the object is to be loaded. The function does not check if the loaded object is of the same name.
#' @param fileName A filename from which the object should be loaded. Usually omitted and constructed automatically.
loadObj<-function(x,fileName=NULL){
  objName <- deparse(substitute(x))
  if(is.null(fileName)&(!is.null(objName))) fileName<-paste(objName,".Rdata",sep="")
  load(file = fileName,.GlobalEnv)
}

#' For a vector x, it computes x[1]/x.
#'
#' @param x A numeric vector. Should not contain 0s.
#' @return A vector computed as x[1]/x.
relInv<-function(x)x[1]/x


#' For a vector x, it computes x[1]/x. If result is not finite (e.g. if certain elements of x are 0), the result is 0.
#'
#' @param x A numeric vector. Can contain 0s.
#' @return A vector computed as x[1]/x. If result is not finite (e.g. if certain elements of x are 0), the result is 0.
relInv2<-function(x){
  x<-x[1]/x
  x[!is.finite(x)]<-0
  x
}

#' Computes sum of squared deviations from the mean using only valid (non NA) values.
#'
#' @param x A vector (or something that can be coerced to a vector) of values (which can include NAs)
#' @return Sum of squared deviations from the mean using only valid (non NA) values.
ssNa<-function(x)ss(na.omit(as.vector(x)))


#' Expands a square matrix by repeating each row/column the specified number of times
#'
#' @param mat A square matrix to be exapanded
#' @param nn A vector of number of times each row/column must be repeated. Its length must match the number of rows/columns
#' @return Sum of squared deviations from the mean using only valid (non NA) values.
expandMat<-function(mat, nn){
    v<-rep(1:length(nn),nn)
  mat[v,v]
}

#' Performed first k-means approach folowed by a homogeneity sum of squared generalized blockmodeling approach using only complete blocks
#'
#' @inheritParams kmBlockORP
#' @param orFun kmBlock currently requires M to be a matrix. If M is a 3 dimmensional array with the third dimmension repesenting relation, the function specified by orFun (must be a function object) is used to combine different relations
#' @param blocks Allowed blocks or pre-specified block image for sum of squares blockmodeling approach. Anyting else but the default probably does not make sense. 
#' @return Sum of squared deviations from the mean using only valid (non NA) values.
kMeansAndSSMlevOR<-function(M, k, n=ncol(M), weightMat=matrix(1), rep=1000, nCores=0, orFun=max, blocks="com", ...){
  
  if(length(dim(M))>2)if(length(dim(M))==3){
    M<-apply(M, c(1,2),sum)
  }else stop("M must have 2 or 3 dimmensions")
  tmpW<-expandMat(weightMat,n)
  
  posWeights<-expandMat(weightMat,k)
  
  kmRes<-kmBlockORP(M = M,k = k,n = n,weights = tmpW, rep = rep, nCores = nCores, ...)
  tclu<-splitCluRes(kmRes)
  if(length(tclu)==1) tclu<-unlist(tclu)
  ssRes<-optParC(M = M, clu=tclu, approach="hom",blocks=blocks,posWeights = posWeights,useMulti = TRUE)
  ssRes$kmRes<-kmRes
  return(ssRes)
}



#' Replaces NaN values by the speficied values (0 by default)
#'
#' @param x A vector or similar where the NaNs are to be replaced.
#' @param rep A value that should replace the NaNs (0 by default).
#' @return x with NaNs replaced by rep.
nanRep<-function(x, rep=0){
  x[is.nan(x)]<-rep
  return(x)
}



#' Computes weights for parts of the multilevel network based on random errors using the SS approach with complete blocks only (compatible with k-means)
#'
#' @param mlNet A multilevel/linked network
#' @param cluParts A partition spliting the units into different sets
#' @param k A vecotor of number of clusters for each set of units in the network.
#' @param mWeights The number of repetitions for computing random errors. Defaults to 1000
#' @return Weights and "intermediate results":
#'  \item{"errArr"}{A 3d array of errors (\code{mWeights} for each part of the network)}
#'  \item{"errMatSum"}{\code{errArr} summed over all repretitions.}
#'  \item{"weightsMat"}{A matrix of weights, one for each part. An inverse of \code{errMatSum} with NaNs replaced by zeros.}
weightsMlSS<-function(mlNet,cluParts, k, mWeights=1000){
  cluParts<-as.numeric(factor(cluParts))
  nn<-table(cluParts)
  mlOrNet<-apply(mlNet, c(1,2),sum)
  parts<-fun.by.blocks(mlOrNet, clu = cluParts, ignore.diag = FALSE, FUN = ss)>0
  errArr<-array(NA,dim=c(dim(parts),mWeights))
  
  for(i1 in 1:dim(parts)[1]){
    for(i2 in 1:dim(parts)[2]){
      if(is.finite(parts[i1,i2])&parts[i1,i2]){
        net<-mlOrNet[cluParts==i1, cluParts==i2]
        if(i1!=i2){
          tnet<-matrix(0, nrow=(nn[i1] + nn[i2]), ncol=(nn[i1] + nn[i2]))
          tnet[cluParts[cluParts %in% c(i1,i2)]==i1, cluParts[cluParts %in% c(i1,i2)]==i2]<-net
          net<-tnet
        }
        for(i in 1:mWeights){
          if(i1==i2){
            tclu<-genRandomPar(k = k[i1],n = nn[i1],addParam=list(genPajekPar = FALSE))
          } else {
            tclu<-genRandomPar(k = k[c(i1,i2)],n = nn[c(i1,i2)],addParam=list(genPajekPar = FALSE))
            tclu[[2]]<-tclu[[2]]+k[i1]
            tclu<-unlist(tclu)
          }
          tCF<-critFunC(net, clu=tclu, approach="hom",homFun="ss", blocks="com")
          
          errArr[i1,i2,i]<-tCF$err
        }
      }
      cat("i1 = ", i1, ", i2 = ", i2,"\n")
    }
  }
  errMatSum<-apply(errArr,c(1,2),sum)
  weightsMat<-relInv2(errMatSum)
  res<-list(errArr=errArr, errMatSum=errMatSum, weightsMat=weightsMat)
}





#' Computes errors for parts of the multilevel network using the SS approach with complete blocks only (compatible with k-means)
#'
#' @param mlNet A multilevel/linked network
#' @param cluParts A partition spliting the units into different sets.
#' @param nn A vector of mumber of units in each set. Either this or \code{cluParts} must be suplied
#' @param clu A partition for which the errors by parts should be comptued.
#' @return \code{x} with NaNs replaced by rep.
errByPartsMlSS<-function(mlNet,cluParts=NULL, nn=NULL, clu){
  if(is.null(cluParts)){
    if(is.null(nn)) stop("Either cluParts or nn must be supplied")
    cluParts<-rep(1:length(nn),times=nn)
  }
  if(!is.list(clu)){
    clu<-split(clu,cluParts)
  }
  cluParts<-as.numeric(factor(cluParts))
  nn<-table(cluParts)
  mlOrNet<-apply(mlNet, c(1,2),sum)
  parts<-fun.by.blocks(mlOrNet, clu = cluParts, ignore.diag = FALSE, FUN = ss)>0
  errMat<-parts
  errMat[,]<-NA
  
  for(i1 in 1:dim(parts)[1]){
    for(i2 in 1:dim(parts)[2]){
      if(is.finite(parts[i1,i2])&parts[i1,i2]){
        net<-mlOrNet[cluParts==i1, cluParts==i2]
        if(i1!=i2){
          tnet<-matrix(0, nrow=(nn[i1] + nn[i2]), ncol=(nn[i1] + nn[i2]))
          tnet[cluParts[cluParts %in% c(i1,i2)]==i1, cluParts[cluParts %in% c(i1,i2)]==i2]<-net
          net<-tnet
        }
          if(i1==i2){
          tclu<-clu[[i1]]
          tclu<-as.numeric(factor(tclu))
        } else {
          tclu<-clu[c(i1,i2)]
          
          tclu<-lapply(tclu,function(x)as.numeric(factor(x)))
          tclu[[2]]<-tclu[[2]]+max(tclu[[1]])
          tclu<-unlist(tclu)
        }
        tCF<-critFunC(net, clu=tclu, approach="hom",homFun="ss", blocks="com")
        
        errMat[i1,i2]<-tCF$err
    }
      #cat("i1 = ", i1, ", i2 = ", i2,"\n")
    }
  }
  return(errMat)
}


#' Optimizes the partitions for parts of the multilevel network using the k-means approach folowed by SS approach with complete blocks only (compatible with k-means)
#'
#' @param mlNet A multilevel/linked network
#' @param cluParts A partition spliting the units into different sets.
#' @param nn A vector of mumber of units in each set. Either this or \code{cluParts} must be suplied
#' @param k A vector of the number of clusters (one for each set of units).
#' @param rep The number of random restarts to try.
#' @return \code{x} with NaNs replaced by rep.
optErrByPartsMlSS<-function(mlNet,cluParts=NULL, nn=NULL, k, rep=100, nCores=0){
  if(is.null(cluParts)){
    if(is.null(nn)) stop("Either cluParts or nn must be supplied")
    cluParts<-rep(1:length(nn),times=nn)
  }
  clu<-cluParts
  clu<-as.numeric(factor(clu))
  nn<-table(clu)
  mlOrNet<-apply(mlNet, c(1,2),sum)
  parts<-fun.by.blocks(mlOrNet, clu = clu, ignore.diag = FALSE, FUN = ss)>0
  errMat<-parts
  errMat[,]<-NA
  
  for(i1 in 1:dim(parts)[1]){
    for(i2 in 1:dim(parts)[2]){
      if(is.finite(parts[i1,i2])&parts[i1,i2]){
        net<-mlOrNet[clu==i1, clu==i2]
      
        if(i1!=i2){
          tnet<-matrix(0, nrow=(nn[i1] + nn[i2]), ncol=(nn[i1] + nn[i2]))
          tnet[clu[clu %in% c(i1,i2)]==i1, clu[clu %in% c(i1,i2)]==i2]<-net
          net<-tnet
          
          tmpK<-k[c(i1,i2)]
          
          n = nn[c(i1,i2)]
        } else {
          tmpK<-k[i1]
          n = nn[i1]
        }
        kmRes<-kmBlockORP(M = net,k = tmpK,n = n, rep = rep, nCores = nCores)
        
        if(i1==i2){
          tclu<-clu(kmRes)  
        } else tclu<-splitCluRes(kmRes)
        
        ssRes<-optParC(M = net, clu=tclu, approach="hom", blocks="com", useMulti = TRUE)
        
          errMat[i1,i2]<-ssRes$err
        
      }
      cat("i1 = ", i1, ", i2 = ", i2,"\n")
    }
  }
  return(errMat)}
