"clu" <-
function(res,which=1,...){
  if("clu" %in% names(res)){
    res$clu
  }else res$best[[which]]$clu
}

"partitions" <- 
function(res)lapply(res$best,function(x)x$clu)
