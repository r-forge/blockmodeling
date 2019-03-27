# to do - here and in C-functions --> put functional blocks before regular !!!

cStatus<-list(
    blockTypes=c("nul", "com", "cfn", "rfn", "reg", "cre", "rre", "avg", "dnc"), #add before average 
    regFuns=c("max","sum","mean"), 
    homFuns=c("ss", "ad"), 
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



########## warning -- this functions needs to be corrected to be more similar to optParC and optRandParC
critFunC<-function(M, clu, approaches, blocks, isTwoMode = NULL, isSym = NULL,
                   diag = 1, IM = NULL, EM = NULL, Earr = NULL, justChange = FALSE, 
                   rowCluChange = c(0, 0), colCluChange = c(0, 0), sameIM = FALSE, 
                   regFun = "max", homFun = "ss", usePreSpecM = NULL, preSpecM = NULL, 
                   save.initial.param = TRUE, relWeights = 1, posWeights = 1, 
                   blockTypeWeights = 1, combWeights = NULL, returnEnv = FALSE){
    if(save.initial.param){
        initial.param<-list(initial.param=tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")))   #saves the inital parameters
    }else initial.param<-NULL
    

    if(length(dim(M))==2) M<-array(M,dim=c(dim(M),length(approaches)))
    #M[,,approaches=="bin"]<-(M[,,approaches=="bin"]>0)*1
    dM<-dim(M)
    if(is.null(isTwoMode)) isTwoMode<-is.list(clu)

    if(!is.list(clu))clu<-list(clu,clu)
    orgClu<-clu
    clu<-lapply(clu,function(x)as.integer(as.factor(x)))
    nUnitsInRCclu<-lapply(clu,function(x)as.integer(table(x)))
    nRCclu<-sapply(nUnitsInRCclu,length)

    # if(is.null(nMode)) nMode<-ifelse(is.list(clu),length(clu),1)
    # if(nMode>1){
        # tmNclu<-sapply(clu,max)
        # for(iMode in 2:nMode){
            # clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
        # }
        
        # clu<-unlist(clu)    
    # }
        
    rowParArr<-matrix(as.integer(0),nrow=dM[1],ncol=nRCclu[1])
    for(i in 1:nRCclu[[1]]){
        rowParArr[1:nUnitsInRCclu[[1]][i],i]<-as.integer(which(clu[[1]]==i)-1)
    }
    colParArr<-matrix(as.integer(0),nrow=dM[2],ncol=nRCclu[2])
    for(i in 1:nRCclu[[2]]){
        colParArr[1:nUnitsInRCclu[[2]][i],i]<-as.integer(which(clu[[2]]==i)-1)
    }


    if(is.null(isSym)){
        isSym<-integer(dM[3])
        if(isTwoMode) {
            isSym[]<-FALSE
        } else {
            for(i in 1:dM[3]) isSym[i]<-all(M[,,i]==t(M[,,i]))
        }
    } else if(length(isSym)==1) isSym<-rep(isSym, dM[3])

    if(isTwoMode)diag<-FALSE
    if(length(diag)!=dM[3]) diag<-rep(diag[1], dM[3])
    if(length(approaches)!=dM[3]&&(length(approaches)==1)) approaches<-rep(approaches[1], dM[3])

    if(is.list(blocks)){
        if(length(blocks)!=dM[3]) stop("the number of relations implied by 'blocks' and by 'M' does not match")
        maxBlockTypes<- max(sapply(blocks,length))
        blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
        for(i in 1:dM[3]){
            nBT<-length(blocks[[i]])
            blocksArr[1:nBT,i,,]<-array(blocks[[i]],dim=c(nBT,nRCclu))
        }
        blocks <- blocksArr
    } else if(is.vector(blocks)){
        maxBlockTypes<-length(blocks)
        blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
        blocksArr[1:length(blocks),,,]<-blocks
        blocks <- blocksArr
    } else if(!is.array(blocks)){
        stop("'blocks' argument should be a vector, a list or an array with appropriate dimmensions")
    }else {
        if(length(dim(blocks))==4){
            maxBlockTypes<-dim(blocks)[1]
            if(any(dim(blocks)!=c(maxBlockTypes,dM[3],nRCclu))) stop("array ('blocks' argument) has a wrong dimensions of dimmensions")
        } else if(length(dim(blocks))==3){
            maxBlockTypes<-dim(blocks)[1]
            blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
            for(i in 1:dM[3]){
                blocksArr[,i,,]<-blocks
            }
            blocks <- blocksArr
        } else if(length(dim(blocks))==2){
            maxBlockTypes<-1
            blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
            for(i in 1:dM[3]){
                blocksArr[1,i,,]<-blocks
            }
            blocks <- blocksArr
        } else stop("array ('blocks' argument) has a wrong number of dimmensions")
    }
    dB<-dim(blocks)

    if(dB[2]!=dM[3])stop("the number of relations implied by 'blocks' and by 'M' does not match")
    if(!all(dB[3:4]==nRCclu))stop("number of clusters implied by 'blocks' and by 'clu' does not match")
    nBlockTypeByBlock<-apply(!is.na(blocks),c(2,3,4),sum)
    blocks[blocks=="null"]<-"nul"
	blocks[blocks=="den"]<-"avg"

    if(is.null(IM)){
        IM<-array(as.integer(99),dim=dB[2:4])
    }else if (length(dim(IM))==2){
        IM<-array(as.integer(factor(IM,levels=cStatus$blockTypes))-as.integer(1),dim=c(dM[3],nRCclu))
    }else{
        IM<-array(as.integer(factor(IM,levels=cStatus$blockTypes))-as.integer(1),dim=dim(IM))
    }

    if(is.null(EM)){
        EM<-array(as.double(Inf),dim=dB[2:4])
    } else EM<-array(as.double(EM),dim=dim(EM))
    if(is.null(Earr)){
        Earr<-array(as.double(Inf),dim=dB)
    }else Earr<-array(as.double(Earr),dim=dim(Earr))
    
    if(length(homFun)==1 & dM[3]>1) homFun<-rep(homFun,dM[3])
	
	  homFun[approaches=="ss"]<-"ss"
	  homFun[approaches=="ad"]<-"ad"
	  approaches[approaches%in%c("ss","ad")]<-"hom"
	  
    homFun<-as.integer(factor(homFun,levels=cStatus$homFuns))-as.integer(1)
    
    
    regFun<-as.integer(factor(regFun,levels=cStatus$regFuns))-as.integer(1)
    if(is.vector(regFun)){
        if(length(regFun)==1){
            regFun <- array(as.integer(regFun),dim=dB)
        }else if (dB[2]==1){
            if(length(regFun)==dB[1]){
                regFunArr <- array(as.integer(NA),dim=dB)
                regFunArr[,,,]<-regFun
                regFun<-regFunArr
            } else stop("'regFun' is a vector of unapropriate length")
        } else if(length(regFun)==dB[2]){
                regFunArr <- array(as.integer(NA),dim=dB)
                regFunArr[,,,]<-regFun
                regFun<-regFunArr
        } else stop("'regFun' is a vector of unapropriate length")
    } else if(is.array(regFun)){
        if(dim(regFun)!=dB){
            stop("'regFun' is an array - dimensions of 'regFun' and 'blocks' do not match")
        }
    } else stop("'regFun' is neither a vector or an array")

    preSpecM<-formatPreSpecM(preSpecMorg=preSpecM,dB=dB,blocks=blocks)

    usePreSpecM<-formatUsePreSpecM(usePreSpecMorg=usePreSpecM,preSpecM=preSpecM,dB=dB,blocks=blocks)

    if(any(approaches=="bin") && (!all(M[,,approaches=="bin"] %in% c(0,1)))){
        for(i in 1:length(approaches)){
            if(approaches[i]=="bin"){
                if(!all(M[,,i] %in% c(0,1))){
                    tmpPreSpecM<-preSpecM[,i,,]
                    if(all(is.na(tmpPreSpecM))){
                        M[,,i]<-(M[,,i]>0)*1
                    } else if(all(tmpPreSpecM==tmpPreSpecM[1,1,1])){
                        M[,,i]<-(M[,,i]>=tmpPreSpecM[1,1,1])*1
                    } else stop("Relation ",i," is not binary but suplied to binary blockmodeling without suitable value in 'preSpecM'!",sep="")
                }
            }
        }
    }
    
    approaches <- as.integer(factor(approaches,levels=cStatus$implementedApproaches))-as.integer(1)
    
    
    

    combWeights<-computeCombWeights(combWeights, dB, blocks, relWeights, posWeights, blockTypeWeights)
    blocks<-array(as.integer(factor(blocks,levels=cStatus$blockTypes)),dim=dim(blocks))-as.integer(1)
    
    M<-apply(M,c(2,3),as.double)
    
    resC<-.C("critFun", M=M, nr=dM[1], nc=dM[2], nRel=dM[3], isTwoMode=as.integer(isTwoMode), isSym=as.integer(isSym), diag=as.integer(diag), nColClus=nRCclu[2], nRowClus=nRCclu[1], nUnitsRowClu=nUnitsInRCclu[[1]], nUnitsColClu=nUnitsInRCclu[[2]], rowParArr=rowParArr, colParArr=colParArr, approaches=approaches, maxBlockTypes=as.integer(maxBlockTypes), nBlockTypeByBlock=array(as.integer(nBlockTypeByBlock),dim=dim(nBlockTypeByBlock)), blocks=blocks, IM=IM, EM=EM, Earr=Earr, err=sum(EM), justChange=as.integer(justChange), rowCluChange=as.integer(rowCluChange), colCluChange=as.integer(colCluChange), sameIM=as.integer(sameIM), regFun=regFun, homFun=homFun, usePreSpec=usePreSpecM, preSpecM=preSpecM,combWeights=combWeights,NAOK=TRUE)
    
    
    res<-c(list(M=M), resC[c("err","EM","Earr")], list(IM=IMaddNames(resC$IM)), list(clu=orgClu), initial.param, list(call=match.call()), if(returnEnv)list(env= environment()) else NULL)
    class(res)<-"crit.fun"
    return(res)
}


optParC<-function(M, clu, approaches, blocks, nMode=NULL,isSym=NULL,diag=1, useMulti=FALSE, maxPar=50, IM=NULL,EM=NULL,Earr=NULL, justChange=TRUE, justMove=FALSE, sameIM=FALSE, regFun="max", homFun = "ss", usePreSpecM = NULL, preSpecM=NULL, minUnitsRowCluster = 1, minUnitsColCluster = 1, maxUnitsRowCluster = 9999, maxUnitsColCluster = 9999, relWeights=1, posWeights=1, blockTypeWeights=1,combWeights=NULL, exchageClusters="all",save.initial.param=TRUE){

    if(save.initial.param){
        initial.param<-list(initial.param=tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")))   #saves the inital parameters
    }else initial.param<-NULL
    

    if(length(dim(M))==2) M<-array(M,dim=c(dim(M),length(approaches)))
    dM<-dim(M)
    if(is.null(nMode)) nMode<-ifelse(is.list(clu),length(clu),1)
    
    if(nMode>1){
        tmN<-sapply(clu,length)
        clu<-lapply(clu,function(x)as.integer(as.factor(x)))
        tmNclu<-sapply(clu,max)
        for(iMode in 2:nMode){
            clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
        }
        
        clu<-unlist(clu)
        if(dM[1]!=length(clu)|dM[2]!=length(clu)){
            warning("Two (and more) mode networks implemented through one mode networks!\nOnly partition, network and blocks arguments are converted if needed!\nIf usePrespecVal and similar arguments are arrays they must be in appropriate format - one mode network with two-mode network in upper right quadrant")
            #currently two mode networks are treated as a special case of one mode networks where 3 "quadrants" of the network are filled with zeros
            oldM<-M
            oldDM<-dim(oldM)
            nUnitsTmp<-length(clu)
            M<-array(0,dim=c(nUnitsTmp,nUnitsTmp,length(approaches)))
            M[1:oldDM[1],((oldDM[1]+1):nUnitsTmp),]<-oldM
            dM<-dim(M)
        }
    }
    
    if(!is.list(clu))clu<-list(clu,clu)
    clu<-lapply(clu,function(x)as.integer(as.factor(x))-as.integer(1))
    nUnitsInRCclu<-lapply(clu,function(x)as.integer(table(x)))
    nRCclu<-sapply(nUnitsInRCclu,length)
    rowParArr<-matrix(as.integer(0),nrow=dM[1],ncol=nRCclu[1])
    for(i in 1:nRCclu[1]){
        rowParArr[1:nUnitsInRCclu[[1]][i],i]<-as.integer(which(clu[[1]]==(i-1))-1)
    }
    colParArr<-matrix(as.integer(0),nrow=dM[2],ncol=nRCclu[2])
    for(i in 1:nRCclu[2]){
        colParArr[1:nUnitsInRCclu[[2]][i],i]<-as.integer(which(clu[[2]]==(i-1))-1)
    }

    if(exchageClusters=="all"){
        if(nMode>1){
            exchageClusters=matrix(as.integer(0),nrow=nRCclu[1],ncol=nRCclu[2])
            
            tmp<-c(0,tmNclu)
            for(imodeNclu in seq_along(tmNclu)){
                tmpInd<-(sum(tmp[1:imodeNclu])+1):sum(tmNclu[1:imodeNclu])
                exchageClusters[tmpInd,tmpInd]=as.integer(1)
            }
        } else{
            exchageClusters=matrix(as.integer(1),nrow=nRCclu[1],ncol=nRCclu[2])
        }
    }
    
    if(is.null(isSym)){
        isSym<-integer(dM[3])
        for(i in 1:dM[3]) isSym[i]<-all(M[,,i]==t(M[,,i]))
    } else if(length(isSym)==1) isSym<-rep(isSym, dM[3])

    #if(isTwoMode)diag<-FALSE #not needed as two mode netowrks are implemented through one-mode networks
    if(length(diag)!=dM[3]) diag<-rep(diag[1], dM[3])
    if(length(approaches)!=dM[3]&&(length(approaches)==1)) approaches<-rep(approaches[1], dM[3])

    if(is.list(blocks)){
        if(length(blocks)!=dM[3]) stop("the number of relations implied by 'blocks' and by 'M' does not match")
        maxBlockTypes<- max(sapply(blocks,length))
        blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
        for(i in 1:dM[3]){
            nBT<-length(blocks[[i]])
            blocksArr[1:nBT,i,,]<-array(blocks[[i]],dim=c(nBT,nRCclu))
        }
        blocks <- blocksArr
    } else if(is.vector(blocks)){
        maxBlockTypes<-length(blocks)
        blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
        blocksArr[1:length(blocks),,,]<-blocks
        blocks <- blocksArr
    } else if(!is.array(blocks)){
        stop("'blocks' argument should be a vector, a list or an array with appropriate dimmensions")
    }else {
        if(length(dim(blocks))==4){
            maxBlockTypes<-dim(blocks)[1]
            if(any(dim(blocks)!=c(maxBlockTypes,dM[3],nRCclu))){
                if(nMode==2){
                    oldBlocks<-blocks
                    blocks<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
                    blocks[,,1:tmNclu[1],(tmNclu[1]+1):sum(tmNclu)]<-oldBlocks
                    blocks[1,,(tmNclu[1]+1):sum(tmNclu),]<-"dnc"
                    blocks[1,,1:tmNclu[1],1:tmNclu[1]]<-"dnc"
                    if(any(dim(blocks)!=c(maxBlockTypes,dM[3],nRCclu))) stop("array ('blocks' argument) has a wrong dimensions of dimensions")    
                } else stop("array ('blocks' argument) has a wrong dimensions of dimensions")
            }
        } else if(length(dim(blocks))==3){
            maxBlockTypes<-dim(blocks)[1]
            blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))      
            if(nMode==2){
                for(i in 1:dM[3]){
                    blocksArr[,i,1:tmNclu[1],(tmNclu[1]+1):sum(tmNclu)]<-blocks 
                }
            } else {
                for(i in 1:dM[3]){
                    blocksArr[,i,,]<-blocks 
                }
            }
            blocks <- blocksArr
            if(nMode==2){
                blocks[1,,(tmNclu[1]+1):sum(tmNclu),]<-"dnc"
                blocks[1,,1:tmNclu[1],1:tmNclu[1]]<-"dnc"
            }
            
        } else if(length(dim(blocks))==2){
            maxBlockTypes<-1
            blocksArr<-array(NA,dim=c(maxBlockTypes,dM[3],nRCclu))
            if(nMode==2){
                for(i in 1:dM[3]){
                    blocksArr[1,i,1:tmNclu[1],(tmNclu[1]+1):sum(tmNclu)]<-blocks
                }
            }else {
                for(i in 1:dM[3]){
                    blocksArr[1,i,,]<-blocks
                }
            }
            blocks<-blocksArr
            if(nMode==2){
                blocks[1,,(tmNclu[1]+1):sum(tmNclu),]<-"dnc"
                blocks[1,,1:tmNclu[1],1:tmNclu[1]]<-"dnc"
            }
        } else stop("array ('blocks' argument) has a wrong number of dimmensions")
    }
    
    dB<-dim(blocks)

    if(dB[2]!=dM[3])stop("the number of relations implied by 'blocks' and by 'M' does not match")
    if(!all(dB[3:4]==nRCclu))stop("number of clusters implied by 'blocks' and by 'clu' does not match")
    nBlockTypeByBlock<-apply(!is.na(blocks),c(2,3,4),sum)
    blocks[blocks=="null"]<-"nul"
	blocks[blocks=="den"]<-"avg"

    if(is.null(IM)){
        IM<-array(as.integer(99),dim=dB[2:4])
    }else if (length(dim(IM))==2){
        IM<-array(as.integer(factor(IM,levels=cStatus$blockTypes))-as.integer(1),dim=c(dM[3],nRCclu))
    }else{
        IM<-array(as.integer(factor(IM,levels=cStatus$blockTypes))-as.integer(1),dim=dim(IM))
    }

    if(is.null(EM)){
        EM<-array(as.double(Inf),dim=dB[2:4])
    } else EM<-array(as.double(EM),dim=dim(EM))
    if(is.null(Earr)){
        Earr<-array(as.double(Inf),dim=dB)
    }else Earr<-array(as.double(Earr),dim=dim(Earr))
    
    if(length(homFun)==1 & dM[3]>1) homFun<-rep(homFun,dM[3])
	
  	homFun[approaches=="ss"]<-"ss"
  	homFun[approaches=="ad"]<-"ad"
  	approaches[approaches%in%c("ss","ad")]<-"hom"
  	
    homFun<-as.integer(factor(homFun,levels=cStatus$homFuns))-as.integer(1)
    
    regFun<-as.integer(factor(regFun,levels=cStatus$regFuns))-as.integer(1)
    if(is.vector(regFun)){
        if(length(regFun)==1){
            regFun <- array(as.integer(regFun),dim=dB)
        }else if (dB[2]==1){
            if(length(regFun)==dB[1]){
                regFunArr <- array(as.integer(NA),dim=dB)
                regFunArr[,,,]<-regFun
                regFun<-regFunArr
            } else stop("'regFun' is a vector of unapropriate length")
        } else if(length(regFun)==dB[2]){
                regFunArr <- array(as.integer(NA),dim=dB)
                regFunArr[,,,]<-regFun
                regFun<-regFunArr
        } else stop("'regFun' is a vector of unapropriate length")
    } else if(is.array(regFun)){
        if(dim(regFun)!=dB){
            stop("'regFun' is an array - dimensions of 'regFun' and 'blocks' do not match")
        }
    } else stop("'regFun' is neither a vector or an array")

    preSpecM<-formatPreSpecM(preSpecMorg=preSpecM,dB=dB,blocks=blocks)
    usePreSpecM<-formatUsePreSpecM(usePreSpecMorg=usePreSpecM,preSpecM=preSpecM,dB=dB,blocks=blocks)

    if(any(approaches=="bin") && (!all(M[,,approaches=="bin"] %in% c(0,1)))){
        for(i in 1:length(approaches)){
            if(approaches[i]=="bin"){
                if(!all(M[,,i] %in% c(0,1))){
                    tmpPreSpecM<-preSpecM[,i,,]
                    if(all(is.na(tmpPreSpecM))){
                        M[,,i]<-(M[,,i]>0)*1
                    } else if(all(tmpPreSpecM==tmpPreSpecM[1,1,1])){
                        M[,,i]<-(M[,,i]>=tmpPreSpecM[1,1,1])*1
                    } else stop("Relation ",i," is not binary but suplied to binary blockmodeling without suitable value in 'preSpec'!",sep="")
                }
            }
        }
    }
    approaches <- as.integer(factor(approaches,levels=cStatus$implementedApproaches))-as.integer(1)
    
    M<-apply(M,c(2,3),as.double)
    
    combWeights<-computeCombWeights(combWeights, dB, blocks, relWeights, posWeights, blockTypeWeights)
	blocks<-array(as.integer(factor(blocks,levels=cStatus$blockTypes)),dim=dim(blocks))-as.integer(1)

																	 
																	 
    
	if(useMulti){		
		bestColParMatrix <- matrix(as.integer(NA),ncol=maxPar,nrow=dM[2])
		bestRowParMatrix <- matrix(as.integer(NA),ncol=maxPar,nrow=dM[1])
		
		resC<-.C("optParMulti", M=M, nr=dM[1], nc=dM[2], nRel=dM[3], isTwoMode= 0 #as.integer(isTwoMode) - two mode networks are currently implemented through onemode networks
		, isSym=as.integer(isSym), diag=as.integer(diag), nColClus=nRCclu[2], nRowClus=nRCclu[1], nUnitsRowClu=nUnitsInRCclu[[1]], nUnitsColClu=nUnitsInRCclu[[2]], rowPar=clu[[1]], colPar=clu[[2]], rowParArr=rowParArr, colParArr=colParArr, approaches=approaches, maxBlockTypes=as.integer(maxBlockTypes), nBlockTypeByBlock=array(as.integer(nBlockTypeByBlock),dim=dim(nBlockTypeByBlock)), blocks=blocks, IM=IM, EM=EM, Earr=Earr, err=sum(EM), justChange=as.integer(justChange), rowCluChange=integer(2), colCluChange=integer(2), sameIM=as.integer(sameIM), regFun=regFun, homFun=homFun, usePreSpec=usePreSpecM, preSpecM=preSpecM, minUnitsRowCluster = as.integer(minUnitsRowCluster), minUnitsColCluster = as.integer(minUnitsColCluster), maxUnitsRowCluster = as.integer(maxUnitsRowCluster), maxUnitsColCluster = as.integer(maxUnitsColCluster), sameErr=as.integer(0), nIter=as.integer(0),combWeights=combWeights, exchageClusters=exchageClusters, maxPar=as.integer(maxPar), bestColParMatrix=bestColParMatrix, bestRowParMatrix=bestRowParMatrix,justMove=justMove, NAOK=TRUE)
        clu<- resC$rowPar
    
	} else{	
		resC<-.C("optPar", M=M, nr=dM[1], nc=dM[2], nRel=dM[3], isTwoMode= 0 #as.integer(isTwoMode) - two mode networks are currently implemented through onemode networks
		, isSym=as.integer(isSym), diag=as.integer(diag), nColClus=nRCclu[2], nRowClus=nRCclu[1], nUnitsRowClu=nUnitsInRCclu[[1]], nUnitsColClu=nUnitsInRCclu[[2]], rowParArr=rowParArr, colParArr=colParArr, approaches=approaches, maxBlockTypes=as.integer(maxBlockTypes), nBlockTypeByBlock=array(as.integer(nBlockTypeByBlock),dim=dim(nBlockTypeByBlock)), blocks=blocks, IM=IM, EM=EM, Earr=Earr, err=sum(EM), justChange=as.integer(justChange), rowCluChange=integer(2), colCluChange=integer(2), sameIM=as.integer(sameIM), regFun=regFun, homFun=homFun, usePreSpec=usePreSpecM, preSpecM=preSpecM, minUnitsRowCluster = as.integer(minUnitsRowCluster), minUnitsColCluster = as.integer(minUnitsColCluster), maxUnitsRowCluster = as.integer(maxUnitsRowCluster), maxUnitsColCluster = as.integer(maxUnitsColCluster), sameErr=as.integer(0), nIter=as.integer(0),combWeights=combWeights,exchageClusters=exchageClusters,justMove=justMove, NAOK=TRUE)
		
		clu<- parArrOne2clu(nUnitsClu=resC$nUnitsRowClu, parArr=resC$rowParArr, nClus=resC$nRowClus)
	}
    
    
#    if(isTwoMode){ # not needed as two-mode networks are implementer through onemode networks
#        clu<- list(
#            parArrOne2clu(nUnitsClu=resC$nUnitsRowClu, parArr=resC$rowParArr, nClus=resC$nRowClus), 
#            parArrOne2clu(nUnitsClu=resC$nUnitsColClu, parArr=resC$colParArr, nClus=resC$nColClus)
#        )
#    } else {
# This (under else) is moved up in to the if(useMulti), as it differs for both functions optPar C functions. Most likely, the below code could be used for both, but is not tested.
#        clu<- parArrOne2clu(nUnitsClu=resC$nUnitsRowClu, parArr=resC$rowParArr, nClus=resC$nRowClus)
#    }

# this is new and experimental
    if(nMode>1){
 		clu<-split(clu, f = rep(1:length(tmN),times=tmN))
		clu<-lapply(clu,function(x)as.integer(as.factor(x)))
        tmNclu<-sapply(clu,max)
        for(iMode in 2:nMode){
            clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
        }
	} else clu<-as.integer(as.factor(clu))
													   
													   
													   
    res<-c(list(M=M), resC[c("err","EM","Earr","sameErr")], list(IM=IMaddNames(resC$IM)), clu=list(clu), initial.param, list(call=match.call()),if(useMulti)list(bestRowParMatrix=bestRowParMatrix),list(resC=resC))
    class(res)<-"optPar"
    return(res)
}











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
return.err=TRUE,#if 'FALSE', only the resoults of crit.fun are returned (a list of all (best) soulutions including errors), else the resoult is list
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
... #paramters to optParC
){
  dots<-list(...) #this might not be need - can be removed and all latter occurencies given sufficent testing. Left for now as there is not enought time.
  if(is.null(switch.names)){
    switch.names<-is.null(blocks)
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
    class(res)<-"opt.more.par"
    return(res)
    })
  
  if(nCores==1||!requireNamespace("doParallel")||!requireNamespace("doRNG")){
      if(nCores!=1) {
        oldWarn<-options("warn")
        options(warn=1)
        warning("Only single core is used as package 'doParallel' or 'doRNG' (or both) is/are not available")
        options(oldWarn)
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
     requireNamespace("doParallel")
     requireNamespace("doRNG")
        if(!getDoParRegistered()){
            if(nCores==0){
                nCores<-detectCores()-1                    
            }
            registerDoParallel(nCores)
        }
        nC<-getDoParWorkers()
        oneRep<-function(i,M,approaches, blocks, n,k,mingr,maxgr,addParam,rep,nC,...){
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
		pkgName<-utils::packageName()
        res<-foreach::foreach(i=1:rep,.combine=c, .packages=pkgName) %dorng% oneRep(i=i,M=M,approaches=approaches, blocks=blocks ,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep,nC=nC,...)
        err<-sapply(res,function(x)x$err)
        nIter<-sapply(res,function(x)x$resC$nIter)
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




