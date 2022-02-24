#' @rdname optRandomParC
#' @param x The result of \code{\link{optRandomParC}}.
#' @method print optMorePar
#' @export
"print.optMorePar" <-
  function(x,...){
  rclu<-clu(x)
  if(is.list(rclu)){
    n<-sapply(rclu,length)
    k<-sapply(rclu,function(x)length(unique(x)))
  } else{
    n<-length(rclu)
    k<-length(unique(rclu))
  } 
  cat("Network size:",sum(n),"\n")
  if(length(n)>1) cat("Network size by sets:",n,"\n")

  if(!is.null(x$initial.param$approaches)){
    cat("\nApproachs (paramter): ")
    if(!is.null(x$initial.param$dots.homFun)){
      cat(paste(x$initial.param$approaches, x$initial.param$dots.homFun, sep = "-"),sep = ", ")  
    } else cat(x$initial.param$approaches,sep = ", ")
  }
  
  if(!is.null(x$initial.param$blocks)){
    cat("\nBlocks (paramter)\n")
    printBlocks(x$initial.param$blocks)
    haveBlocks<-TRUE
  } else haveBlocks<-FALSE
  
  cat("\nSizes of clusters:")
  if(length(n)==1) {
    print(table(clu(x)))
  }else{
    for(i in 1:length(n)){
      cat("Set",i,"\n")
      print(table(clu(x)[[i]]))
    }
  }  
  
  rIM<-IM(x)
  if(haveBlocks){
    printIM<-length(x$initial.param$blocks)>1
  }else{
    printIM<-!all(rIM==rIM[1])
  }
  if(printIM){
    cat("\nIM\n")
    if(length(dim(rIM))>2){
      for(i in 1:dim(rIM)[1]){
        cat("Relation ",i,"\n")
        print(data.frame(rIM[i,,],check.names = FALSE))
      }
    } else print(data.frame(rIM,check.names = FALSE))
  }
  cat("\nError:",min(x$err),"\n")
  if(length(x$best)>1) cat(length(x$best),"solutions with minimal error exits. Only results for the first one are shown above!\n")
}