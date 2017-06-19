"err" <-
function(res,...){
  if(is.null(res[["best"]])){
    min(res$err)
  }else res$best[[1]]$err
}

