cStatus<-list(
    blockTypes=c("nul", "com", "cdo", "rdo", "cfn", "rfn", "reg", "cre", "rre", "avg", "dnc"), 
    regFuns=c("max","sum","mean"), 
    homFuns=c("ss", "ad", "bll"), 
    implementedApproaches=c("hom", "bin","val")
#   ,maxBlockTypes=as.integer(10)
)
# zgornje spremenljivke morajo biti enake kot v C-ju (blockmodeling.c)

allInDimEqual<-function(arr,d)all(apply(arr,d,function(x){x<-as.vector(x);all(x==x[1])}))

clu2parArr<-function(clu){
    if(!is.list(clu))clu<-list(clu,clu)
    nrc<-sapply(clu,length)
    clu<-lapply(clu,function(x)as.integer(as.factor(x)))
    nUnitsInRCclu<-lapply(clu,function(x)as.integer(table(x)))
    nRCclu<-sapply(nUnitsInRCclu,length)
    rowParArr<-matrix(as.integer(0),nrow=nrc[1],ncol=nRCclu[1])
    for(i in clu[[1]]){
        rowParArr[1:nUnitsInRCclu[[1]][i],i]<-as.integer(which(clu[[1]]==i)-1)
    }
    colParArr<-matrix(as.integer(0),nrow=nrc[2],ncol=nRCclu[2])
    for(i in clu[[2]]){
        colParArr[1:nUnitsInRCclu[[2]][i],i]<-as.integer(which(clu[[2]]==i)-1)
    }
    return(list(rowParArr=rowParArr,colParArr=colParArr,nUnitsInRCclu=nUnitsInRCclu, nRCclu=nRCclu, nrc=nrc))
}


parArr2clu<-function(nUnitsRowClu, nUnitsColClu, rowParArr, colParArr, nColClus=NULL, nRowClus=NULL){
    clu<-list(parArrOne2clu(nUnitsClu=nUnitsRowClu, parArr=rowParArr, nClus=nRowClus),parArrOne2clu(nUnitsClu=nUnitsColClu, parArr=colParArr, nClus=nColClus))
}


parArrOne2clu<-function(nUnitsClu, parArr,nClus=NULL){
    if(is.null(nClus)){
        nClus<-dim(parArr)[2]
    } else {
        if(nClus!=dim(parArr)[2]) warning("Number of clusters and dimmension of the partition array do not match")
    }
    n<-sum(nUnitsClu)
    clu<-rep(NA,n)
    for(i in 1:nClus){
        clu[parArr[(1:nUnitsClu[i]),i]+1]<-i
    }
    return(clu)
}


IMaddNames<-function(IM){
    array(factor(IM+1,labels=cStatus$blockTypes,levels=1:length(cStatus$blockTypes)),dim=dim(IM))
}



formatPreSpecM<-function(preSpecMorg,dB,blocks){
    if(is.null(preSpecMorg)){
        preSpecM <- array(as.double(NA),dim=dB)
    } else if (is.vector(preSpecMorg)){
        if(length(preSpecMorg)==1){
			preSpecM <- array(as.double(preSpecMorg),dim=dB)
		} else if(length(preSpecMorg)==dB[2]){
            preSpecM <- array(as.double(NA),dim=dB)
            for(i in 1:dB[2]){
                preSpecM[,i,,]<-as.double(preSpecMorg[i])
            }
        } else if((dB[2]==1) & (length(preSpecMorg)==dB[1]) &  allInDimEqual(blocks,1)){
            preSpecM <- array(as.double(NA),dim=dB)
            for(i in 1:dB[1]){
                preSpecM[i,,,]<-as.double(preSpecMorg[i])
            }           
        } else stop("'",deparse(substitute(preSpecMorg)),"' is a vector with unexpected length")
    } else if(is.array(preSpecMorg)){
        preSpecM <- array(as.double(preSpecMorg),dim=dim(preSpecMorg))
        if(any(dim(preSpecM)!=dB)){
            stop("dimensions of '",deparse(substitute(preSpecMorg)),"' and 'blocks' do not match")
        }   
    }
    return(preSpecM)
}


computeCombWeights<-function(combWeights, dB, blocks, relWeights, posWeights, blockTypeWeights){
    if(!is.null(combWeights)){
        if(all(dim(combWeights)==dB)){
            combWeights<-array(as.double(combWeights),dim=dim(combWeights))
            return(combWeights)
        }
        warning("Dimmensions of the combWeights does not match the dimmensions of blocks!\nIt will not be used!\nIf possible it will be computed using other weights!")
    }
    combWeights<-array(as.double(1),dim=dB)
    
    relWeights<-as.double(relWeights)
    if(length(relWeights)!=dB[2]){
        if(length(relWeights)==1) relWeights<-rep(relWeights,dB[2]) else stop("To relWeights should have length equal to the number of relations!")  
    }
    for(i in 1:dB[2]){
        combWeights[,i,,]<-combWeights[,i,,]*relWeights[i]
    }
    if(all(dim(posWeights)!=dB[3:4])){
        if(length(posWeights)==1) posWeights<-array(posWeights,dim=dB[3:4]) else stop("To posWeights should have the same dimensions as block image!")  
    }
    posWeights<-array(as.double(posWeights), dim=dim(posWeights))
    
    for(i in 1:dB[3]){
        for(j in 1:dB[4]){
            combWeights[,,i,j]<-combWeights[,,i,j]*posWeights[i,j]
        }
    }
    
    if(!(is.numeric(blockTypeWeights)&all(names(blockTypeWeights)%in%cStatus$blockTypes))) stop("blockTypeWeights must be a numeric named vector with names from: ", paste(cStatus$blockTypes, collapse=", "))
    
    for(i in names(blockTypeWeights)){
		tWhich <- blocks==i
		tWhich[is.na(tWhich)]<-FALSE
        combWeights[tWhich]<-blockTypeWeights[i]* combWeights[tWhich]
    }
    return(combWeights)
}


formatUsePreSpecM<-function(usePreSpecMorg,preSpecM,dB,blocks){
    if(is.null(usePreSpecMorg)){
        usePreSpecM<- !is.na(preSpecM)
    }else if(is.vector(usePreSpecMorg)){
        if(length(usePreSpecMorg)==dB[2]){
            usePreSpecM <- array(as.integer(NA),dim=dB)
            for(i in 1:dB[2]){
                usePreSpecM[,i,,]<-as.integer(usePreSpecMorg[i])
            }
        } else if((dB[2]==1) & (length(usePreSpecMorg)==dB[1]) &  allInDimEqual(blocks,1)){
            usePreSpecM <- array(as.integer(NA),dim=dB)
            for(i in 1:dB[1]){
                usePreSpecM[i,,,]<-as.integer(usePreSpecMorg[i])
            }           
        } else stop("'",deparse(substitute(usePreSpecM)),"' is a vector with unexpected length")
    } else if(is.array(usePreSpecMorg)){
        if(any(dim(usePreSpecMorg)!=dB)){
            stop("dimensions of '",deparse(substitute(usePreSpecM)),"' and 'blocks' do not match")
        }   
        usePreSpecM <- array(as.integer(usePreSpecMorg),dim=dim(usePreSpecMorg))
    }
    return(usePreSpecM)
}

