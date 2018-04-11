"gplot1" <-function(M,diag=TRUE,displaylabels=TRUE,boxed.labels=FALSE,loop.cex=4,edge.lwd=1,edge.col="default",rel.thresh=0.05,...){
  if(!requireNamespace("sna", quietly = TRUE)){
    M[M<(max(M)*rel.thresh)]<-0  
    if(edge.col[1]=="default") edge.col<-gray(1-M/max(M))
    edge.col<-edge.col[edge.col!=gray(1)]
    sna::gplot(dat=M,diag=diag,displaylabels=displaylabels,boxed.labels=boxed.labels,loop.cex=loop.cex,edge.lwd=edge.lwd,edge.col=edge.col,...)
  } else stop("Package \"sna\" is needed for this function to work. Please install it.",
                   call. = FALSE)
}


"gplot2" <-
function(M,uselen=TRUE,usecurve=TRUE,edge.len=0.001,diag=TRUE,displaylabels=TRUE,boxed.labels=FALSE,loop.cex=4,arrowhead.cex=2.5,edge.lwd=1,edge.col="default",rel.thresh=0.05,...){
  if(!requireNamespace("sna", quietly = TRUE)){
    M[M<(max(M)*rel.thresh)]<-0
    if(edge.col[1]=="default") edge.col<-gray(1-M/max(M))
    edge.col<-edge.col[edge.col!=gray(1)]
    sna::gplot(dat=M,uselen=uselen,usecurve=usecurve,edge.len=edge.len,diag=diag,displaylabels=displaylabels,boxed.labels=boxed.labels,loop.cex=loop.cex,arrowhead.cex=arrowhead.cex,edge.lwd=edge.lwd,edge.col=edge.col,...)
  } else stop("Package \"sna\" is needed for this function to work. Please install it.",
              call. = FALSE)
}

