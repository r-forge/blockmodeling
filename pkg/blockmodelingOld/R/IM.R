"IM" <-
function(res,which=1,...){
  if(class(res)=="opt.more.par"){
    return(res$best[[which]]$IM)
  } else return(res$IM)
}

