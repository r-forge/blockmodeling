#' @encoding UTF-8
#' @title Functions for Generalized blockmodeling for valued networks
#' 
#' @description
#' Functions for implementation of Generalized blockmodeling for valued
#' networks where the values of the ties are assumed to be measured on at least interval
#' scale. \code{critFunC} calculates the criterion function, based on the network, partition and blockmodel/equivalece.
#' \code{optParC} optimizes a partition based on the criterion function based on a local search algorithm.
#' 
#' @param M A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation. The network can have one or more modes (diferent kinds of units with no ties among themselves). If the network is not two-mode, the matrix must be square.
#' @param clu A partition. Each unique value represents one cluster. If the nework is one-mode, than this should be a vector, else a list of vectors, one for each mode. Similarly, if units are comprised of several sets, \code{clu} should be the list containing one vector for each set.
#' @param approaches One of the approaches (for each relation in multi-relational netowrks in a vector) described in Žiberna (2007). Possible values are:\cr
#' "bin" - binary blockmodeling,\cr
#' "val" - valued blockmodeling,\cr
#' "hom" - homogeneity blockmodeling,\cr
#' "ss" - sum of squares homogeneity blockmodeling, and\cr
#' "ad" - absolute deviations homogeneity blockmodeling.\cr
#' \cr
#' The last two options are "shorthand" for specifying \code{approaches="hom"} and  \code{homFun} to either \code{"ss"} or  \code{"ad"}.
#' @param blocks A vector, a list of vectors or an array with names of allowed blocy types. \cr
#'   \cr
#'   Only listing of allowed block types (blockmodel is not pre-specified).\cr
#'   A vector with names of allowed block types. For multi-relational networks, it can be a list of such vectors. For \code{approaches = "bin"} or \code{approaches = "val"}, at least two should be selected. Possible values are:\cr
#'   \code{"nul"} - null or empty block\cr
#'   \code{"com"} - complete block\cr
#'   \code{"rdo"}, \code{"cdo"} - row and column-dominant blocks (binary and valued approach only)\cr
#'   \code{"reg"} - (f-)regular block\cr
#'   \code{"rre"}, \code{"cre"} - row and column-(f-)regular blocks\cr
#'   \code{"rfn"}, \code{"cfn"} - row and column-dominant blocks (binary, valued only)\cr
#'   \code{"den"} - density block (binary approach only)\cr
#'   \code{"avg"} - average block (valued approach only)\cr
#'   \code{"dnc"} - do not care block - the error is always zero\cr
#'   The ordering is important, since if several block types have identical error, the first on the list is selected.\cr\cr
#'   A pre-specified blockmodel.\cr
#'   An array with four dimensions (see example below). The third and the fourth represent the clusters (for rows and columns). The first is as long as the maximum number of allows block types for a given block. If some block has less possible block types, the empty slots should have values \code{NA}. The second dimension is the number of relations (1 for single-relational networks). The values in the array should be the ones from above. The array can have only three dimensions in case of one-relational networks or if the same pre-specified blockmodel is assumed for all relations. Further, it can have only two dimensions, if in addition only one block type is allowed per block.
#' @param isTwoMode \code{1} for one-mode networks and \code{2} for two-mode networks. The default value is set to \code{NULL}.
#' @param isSym Specifying if the matrix (for each relation) is symmetric.
#' @param diag Should the special status of diagonal be acknowledged. A single number or a vector equal to the number of relation. The default value is set to \code{1}. Codes: \cr
#'    \code{0} - diagonal is treated in the same way as other values \cr
#'    \code{1} - diagonal is  treated separately, or \cr
#'    \code{2} - diagonal values are ignored. \cr
#' @param IM The obtained image for objects. For debugging purposes only.
#' @param EM Block errors by blocks. For debugging purposes only.
#' @param Earr The array of errors for all allowed block types by next dimensions: allowed block types, relations, row clusters and column clusters. The dimensions should match the dimensions of the block argument if specified as an array. For debugging purposes only.
#' @param justChange Value specifying if only the errors for changed clusters should be computed. Used only for debugging purposes by developers.
#' @param rowCluChange An array holding the two row clusters where the change occured. Used only for debugging purposes by developers.
#' @param colCluChange An array holding the col row clusters where the change occured. Used only for debugging purposes by developers.
#' @param sameIM Should we demand the same blockmodel image for all relations. The default value is set to \code{FALSE}.
#' @param regFun Function f used in row-f-regular, column-f-regular, and f-regular blocks. Not used in binary approach. For multi-relational networks, it can be a vector of such character strings. The default value is set to \code{"max"}.
#' @param homFun In case of homogeneity blockmodeling two variability criteria can be used: \code{"ss"} - sum of squares (set by default) and \code{"ad"} -
#' absolute deviations.
#' @param usePreSpecM Specifying weather a pre-specified value should be used when computing inconsistency.
#' @param preSpecM Sufficient value for individual cells for valued approach. Can be a number or a character string giving the name of a function. Set to \code{"max"} for implicit approach. For multi-relational networks, it can be a vector of such values. In case ob binary blockmodeling this argument is a threshold used for binerizing the network. Therefore all values with values lower than \code{preSpecM} are recoded into 0s, all other into 1s. For multi-relational networks, it can be a vector of such values. In case of pre-specified blockmodeling, it can have the same dimensions as \code{blocks}.
#' @param save.initial.param Should the inital parameters (\code{approaches}, ...) be saved. The default value is \code{TRUE}.
#' @param relWeights Weights for all type of relations in a blockmodel. The default value is set to \code{1}.
#' @param posWeights Weigths for positions in the blockmodel (the dimensions must be the same as the error matrix (rows, columns)). For now this is a matix (two-dimensional) even for multi-relational networks.
#' @param blockTypeWeights Weights for each type of block used, if they are to be different across block types (see \code{blocks} above). It must be suplied in form of a named vector, where the names are one or all allowed block types from \code{blocks}. If only some block types are specified, the other have a default weight of 1. The default value is set to \code{1}.
#' @param combWeights Weights for all type of block used, The default value is set to \code{NULL}.The dimension must be the same as \code{blocks}, if \code{blocks} would be specified in array format (which is usual in pre-specified case).
#' @param returnEnv Should the function also return the environment after its completion.
#' @param mulReg Should the errors that apply to rows/columns (and not to cells) should be multiplied by number of rows/columns. Defaults to TRUE.
#' @param addGroupLlErr  Used only when stochastic generalized blockmodeling is used. Should the total error included the part based on sizes of groups. Defaults to TRUE. Will return wrong results for two-mode networks if critFunC is called directly (should be fine if called via optParC function).
#' 
#' 
#'
#' @return 
#' \code{critFunC} returns a list containing:
#' \item{M}{The matrix of the network analyzed.}
#' \item{err}{The error or inconsistency emplirical network with the ideal network for a given blockmodel (model, approach,...) and paritition.}
#' \item{clu}{The analyzed partition.}
#' \item{EM}{Block errors by blocks.}
#' \item{IM}{The obtained image for objects.}
#' \item{BM}{Block means by block - only for Homogeneity blockmodeling.}
#' \item{Earr}{The array of errors for all allowed block types by next dimensions: allowed block types, relations, row clusters and column clusters. The dimensions should match the dimensions of the block argument if specified as an array.}\cr
#' \code{optParC} returns a list containing:
#' \item{M}{The matrix of the network analyzed.}
#' \item{err}{The error or inconsistency emplirical network with the ideal network for a given blockmodel (model, approach,...) and paritition.}
#' \item{clu}{The analyzed partition.}
#' \item{EM}{Block errors by blocks.}
#' \item{IM}{The obtained image for objects.}
#' \item{BM}{Block means by block - only for Homogeneity blockmodeling.}
#' \item{Earr}{The array of errors for all allowed block types by next dimensions: allowed block types, relations, row clusters and column clusters. The dimensions should match the dimensions of the block argument if specified as an array.}
#' \item{useMulti}{The value of the input paramter \code{useMulti}.}
#' \item{bestRowParMatrix}{(If \code{useMulti = TRUE}) Matrix, where there are different solutions for columns, where rows represent units.}
#' \item{sameErr}{The number of partitions with the minimum value of the criterion function.}
#' 
#' @references Doreian, P., Batagelj, V., & Ferligoj, A. (2005). Generalized blockmodeling, (Structural analysis in the social sciences, 25). Cambridge [etc.]: Cambridge University Press.
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2007). Generalized Blockmodeling of Valued Networks. Social Networks, 29(1), 105-126. doi: 10.1016/j.socnet.2006.04.002
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2008). Direct and indirect approaches to blockmodeling of valued networks in terms of regular equivalence. Journal of Mathematical Sociology, 32(1), 57-84. doi: 10.1080/00222500701790207
#' 
#' \enc{Žiberna, A.}{Ziberna, A.} (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002
#'
#' @examples
#' # Generating a simple network corresponding to the simple Sum of squares
#' # Structural equivalence with blockmodel:
#' # nul com
#' # nul nul
#' n <- 20
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(5, 15))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#'
#' # Computation of criterion function with the correct partition
#' res <- critFunC(M = net, clu = clu, approaches = "hom", homFun = "ss", blocks = "com")
#' res$err # The error is relatively small
#' plot(res)
#'
#' # Computation of criterion function with the correct partition and correct pre-specified blockmodel
#' # Prespecified blockmodel used
#' # nul com
#' # nul nul
#' B <- array(NA, dim = c(1, 1, 2, 2))
#' B[1, 1, , ] <- "nul"
#' B[1, 1, 1, 2] <- "com"
#' B[1, 1, , ]
#' res <- critFunC(M = net, clu = clu, approaches = "hom", homFun = "ss", blocks = B)
#' res$err # The error is relatively small
#' res$IM
#' plot(res)
#'
#' # Computation of criterion function with the correct partition
#' # and pre-specified blockmodel with some alternatives
#' # Prespecified blockmodel used
#' # nul nul|com
#' # nul nul
#' B <- array(NA, dim = c(2, 2, 2))
#' B[1, , ] <- "nul"
#' B[2, 1, 2] <- "com"
#' res <- critFunC(M = net, clu = clu, approaches = "hom", homFun = "ss", blocks = B)
#' res$err # The error is relatively small
#' res$IM
#' plot(res)
#'
#' # Computation of criterion function with random partition
#' set.seed(1)
#' clu.rnd <- sample(1:2, size = n, replace = TRUE)
#' res.rnd <- critFunC(M = net, clu = clu.rnd, approaches = "hom",
#' homFun = "ss", blocks = "com")
#' res.rnd$err # The error is larger
#' plot(res.rnd)
#'
#' # Adapt network for Valued blockmodeling with the same model
#' net[net > 4] <- 4
#' net[net < 0] <- 0
#'
#' # Computation of criterion function with the correct partition
#' res <- critFunC(M = net, clu = clu, approaches = "val",
#' blocks = c("nul", "com"), preSpecM = 4)
#' res$err # The error is relatively small
#' res$IM
#' # The image corresponds to the one used for generation of
#' # The network
#' plot(res)
#'
#' # Optimizing one partition
#' res <- optParC(M = net, clu = clu.rnd,
#'    approaches = "hom", homFun = "ss", blocks = "com")
#' plot(res) # Hopefully we get the original partition
#'
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' @seealso \code{\link{optRandomParC}}, \code{\link{IM}}, \code{\link{clu}}, \code{\link{err}}, \code{\link{plot.critFun}}
#' @keywords cluster graphs
#' @import methods
#' 
#' @export

########## warning -- this functions needs to be corrected to be more similar to optParC and optRandParC

critFunC<-function(M, clu, approaches, blocks, isTwoMode = NULL, isSym = NULL,
                   diag = 1, IM = NULL, EM = NULL, Earr = NULL, justChange = FALSE, 
                   rowCluChange = c(0, 0), colCluChange = c(0, 0), sameIM = FALSE, 
                   regFun = "max", homFun = "ss", usePreSpecM = NULL, preSpecM = NULL, 
                   save.initial.param = TRUE, relWeights = 1, posWeights = 1, 
                   blockTypeWeights = 1, combWeights = NULL, returnEnv = FALSE, mulReg=TRUE, addGroupLlErr=TRUE){
  if(save.initial.param){
    initial.param<-list(initial.param=tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")))   #saves the inital parameters
  }else initial.param<-NULL
  
  uniqueBlocks<-unique(unlist(unclass(blocks)))
  if(all(is.na(uniqueBlocks))) stop("No block types are specified!")
  blocksOk<-uniqueBlocks%in%c(blockmodeling:::cStatus$blockTypes,NA)
  if(all(blocksOk)==FALSE) {
	stop("Block types ", paste(uniqueBlocks[!blocksOk], collapse=", ")," are not supported!")
  }
  
  if(length(dim(M))==2) M<-array(M,dim=c(dim(M),length(approaches)))
  #M[,,approaches=="bin"]<-(M[,,approaches=="bin"]>0)*1
  dM<-dim(M)
  if(is.null(isTwoMode)) {
	isTwoMode<-is.list(clu)&& (length(clu)==2)&& (sum(sapply(clu, length))==sum(dim(M)[1:2]))
  }
  
  if(is.list(clu)){
   tmNclu<-sapply(clu,function(x)length(unique(x)))
   tmN<-sapply(clu,length)
  # for(iMode in 2:nMode){
  # clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
  # }
  
  # clu<-unlist(clu)    
  } else {
    tmNclu<-length(unique(clu))
	tmN<-length(clu)
  }  
  
  if(!isTwoMode && is.list(clu)){
	if(sum(sapply(clu, length))==dim(M)[1]) {
		newClu<-c()
		tmpMaxClu<-0
		for(iClu in clu){
			iClu<- as.integer(as.factor(iClu))+tmpMaxClu
			tmpMaxClu <- max(iClu)
			newClu<-c(newClu,iClu)
		}
		clu <- newClu
	} else {
		stop("Clu does not seem to be compatible with M!")
	}
  }

  if(!is.list(clu))clu<-list(clu,clu)
  orgClu<-clu
  clu<-lapply(clu,function(x)as.integer(as.factor(x)))
  nUnitsInRCclu<-lapply(clu,function(x)as.integer(table(x)))
  nRCclu<-sapply(nUnitsInRCclu,length)
  
  if(addGroupLlErr && homFun=="bll"&&(!isTwoMode)){
	nrInSetByClusters<-rep(tmN, tmNclu)
  } else {
	nrInSetByClusters <- as.double(rep(0,nRCclu[1]))
  }
     
  

 
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
          } else if(all(tmpPreSpecM==tmpPreSpecM[1])){
            M[,,i]<-(M[,,i]>=tmpPreSpecM[1])*1
          } else stop("Relation ",i," is not binary but supplied to binary blockmodeling without suitable value in 'preSpecM'!",sep="")
        }
      }
    }
  }
  
  approaches <- as.integer(factor(approaches,levels=cStatus$implementedApproaches))-as.integer(1)
  
  
  
  
  combWeights<-computeCombWeights(combWeights, dB, blocks, relWeights, posWeights, blockTypeWeights)
  blocks<-array(as.integer(factor(blocks,levels=cStatus$blockTypes)),dim=dim(blocks))-as.integer(1)
  
  M<-apply(M,c(2,3),as.double)
  

  
  resC<-.C("critFun", M=M, nr=dM[1], nc=dM[2], nRel=dM[3], isTwoMode=as.integer(isTwoMode), isSym=as.integer(isSym), diag=as.integer(diag), nColClus=nRCclu[2], nRowClus=nRCclu[1], nUnitsRowClu=nUnitsInRCclu[[1]], nUnitsColClu=nUnitsInRCclu[[2]], rowParArr=rowParArr, colParArr=colParArr, approaches=approaches, maxBlockTypes=as.integer(maxBlockTypes), nBlockTypeByBlock=array(as.integer(nBlockTypeByBlock),dim=dim(nBlockTypeByBlock)), blocks=blocks, IM=IM, EM=EM, Earr=Earr, err=sum(EM), justChange=as.integer(justChange), rowCluChange=as.integer(rowCluChange), colCluChange=as.integer(colCluChange), sameIM=as.integer(sameIM), regFun=regFun, homFun=homFun, usePreSpec=usePreSpecM, preSpecM=preSpecM,combWeights=combWeights, mulReg=as.integer(mulReg), nrInSetByClusters= as.integer(nrInSetByClusters), NAOK=TRUE)
  
  
  res<-c(list(M=M), resC[c("err","EM","Earr")], list(IM=IMaddNames(resC$IM)), list(clu=orgClu), initial.param, list(call=match.call()), if(returnEnv)list(env= environment()) else NULL)
  class(res)<-"critFun"
  return(res)
}


 
#' @rdname critFunC
#' 
#' @param nMode Number of nodes. If \code{NULL}, then determined from \code{clu}.
#' @param useMulti Which version of local search should be used. The default value is set to \code{FALSE}. If \code{FALSE}, first possible all moves in random order and then all possible exchanges in random order are tired. When a move with lower value of criterion function is found, the algorithm moves to this new partition. If \code{TRUE} the version of local search where all possible moves and exchanges are tired first and then the one with the lowest error is selected and used. In this case, several optimal partitions are found. \code{maxPar} best partitions are returned. 
#' @param maxPar The number of partitions with optimal criterion fuction to be returned. Only used If \code{useMulti} is \code{TRUE}.
#' @param minUnitsRowCluster Minimum number of units in row cluster.
#' @param minUnitsColCluster Minimum number of units in col cluster.
#' @param maxUnitsRowCluster Maximum number of units in row cluster.
#' @param maxUnitsColCluster Maximum number of units in col cluster.
#' @param exchageClusters A matrix of dimensions "number of clusters" x "number of clusters" indicating to which clusters can units from a specific cluster be moved. Useful for multilevel blockmodeling or/in some other cases where some units cannot mix.
#' 
#' @export


optParC<-function(M, clu, approaches, blocks, nMode=NULL,isSym=NULL,diag=1, useMulti=FALSE, maxPar=50, IM=NULL,EM=NULL,Earr=NULL, justChange=TRUE, sameIM=FALSE, regFun="max", homFun = "ss", usePreSpecM = NULL, preSpecM=NULL, minUnitsRowCluster = 1, minUnitsColCluster = 1, maxUnitsRowCluster = 9999, maxUnitsColCluster = 9999, relWeights=1, posWeights=1, blockTypeWeights=1,combWeights=NULL, exchageClusters="all",save.initial.param=TRUE, mulReg=TRUE, addGroupLlErr=TRUE){
  
  if(save.initial.param){
    initial.param<-list(initial.param=tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")))   #saves the inital parameters
  }else initial.param<-NULL
  
  uniqueBlocks<-unique(unlist(unclass(blocks)))
  if(all(is.na(uniqueBlocks))) stop("No block types are specified!")
  blocksOk<-uniqueBlocks%in%c(blockmodeling:::cStatus$blockTypes,NA)
  if(all(blocksOk)==FALSE) {
	stop("Block types ", paste(uniqueBlocks[!blocksOk], collapse=", ")," are not supported!")
  }
  
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
  } else {
	tmN<-length(clu)
	tmNclu<-length(unique(clu))
  }
  
  if(addGroupLlErr && homFun=="bll"){
	nrInSetByClusters<-rep(tmN, tmNclu)
  } else {
	nrInSetByClusters <- as.double(rep(0,nRCclu[1]))
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
    } else stop("'regFun' is a vector of inappropriate length")
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
          } else if(all(tmpPreSpecM==tmpPreSpecM[1])){
            M[,,i]<-(M[,,i]>=tmpPreSpecM[1])*1
          } else stop("Relation ",i," is not binary but supplied to binary blockmodeling without suitable value in 'preSpec'!",sep="")
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
             , isSym=as.integer(isSym), diag=as.integer(diag), nColClus=nRCclu[2], nRowClus=nRCclu[1], nUnitsRowClu=nUnitsInRCclu[[1]], nUnitsColClu=nUnitsInRCclu[[2]], rowPar=clu[[1]], colPar=clu[[2]], rowParArr=rowParArr, colParArr=colParArr, approaches=approaches, maxBlockTypes=as.integer(maxBlockTypes), nBlockTypeByBlock=array(as.integer(nBlockTypeByBlock),dim=dim(nBlockTypeByBlock)), blocks=blocks, IM=IM, EM=EM, Earr=Earr, err=sum(EM), justChange=as.integer(justChange), rowCluChange=integer(2), colCluChange=integer(2), sameIM=as.integer(sameIM), regFun=regFun, homFun=homFun, usePreSpec=usePreSpecM, preSpecM=preSpecM, minUnitsRowCluster = as.integer(minUnitsRowCluster), minUnitsColCluster = as.integer(minUnitsColCluster), maxUnitsRowCluster = as.integer(maxUnitsRowCluster), maxUnitsColCluster = as.integer(maxUnitsColCluster), sameErr=as.integer(0), nIter=as.integer(0),combWeights=combWeights, exchageClusters=exchageClusters, maxPar=as.integer(maxPar), bestColParMatrix=bestColParMatrix, bestRowParMatrix=bestRowParMatrix, mulReg=as.integer(mulReg), nrInSetByClusters = as.integer(nrInSetByClusters), NAOK=TRUE)
    clu<- resC$rowPar
    
  } else{	
    resC<-.C("optPar", M=M, nr=dM[1], nc=dM[2], nRel=dM[3], isTwoMode= 0 #as.integer(isTwoMode) - two mode networks are currently implemented through onemode networks
             , isSym=as.integer(isSym), diag=as.integer(diag), nColClus=nRCclu[2], nRowClus=nRCclu[1], nUnitsRowClu=nUnitsInRCclu[[1]], nUnitsColClu=nUnitsInRCclu[[2]], rowParArr=rowParArr, colParArr=colParArr, approaches=approaches, maxBlockTypes=as.integer(maxBlockTypes), nBlockTypeByBlock=array(as.integer(nBlockTypeByBlock),dim=dim(nBlockTypeByBlock)), blocks=blocks, IM=IM, EM=EM, Earr=Earr, err=sum(EM), justChange=as.integer(justChange), rowCluChange=integer(2), colCluChange=integer(2), sameIM=as.integer(sameIM), regFun=regFun, homFun=homFun, usePreSpec=usePreSpecM, preSpecM=preSpecM, minUnitsRowCluster = as.integer(minUnitsRowCluster), minUnitsColCluster = as.integer(minUnitsColCluster), maxUnitsRowCluster = as.integer(maxUnitsRowCluster), maxUnitsColCluster = as.integer(maxUnitsColCluster), sameErr=as.integer(0), nIter=as.integer(0),combWeights=combWeights,exchageClusters=exchageClusters, mulReg=as.integer(mulReg), nrInSetByClusters = as.integer(nrInSetByClusters), NAOK=TRUE)
    
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
  } else {
	clu<-as.integer(as.factor(clu))
  }
  
  
  
  res<-c(list(M=M), resC[c("err","EM","Earr","sameErr")], list(IM=IMaddNames(resC$IM)), clu=list(clu), initial.param, list(call=match.call()),if(useMulti)list(bestRowParMatrix=bestRowParMatrix),list(resC=resC))
  class(res)<-"optPar"
  return(res)
}

#' @useDynLib blockmodeling, .registration = TRUE