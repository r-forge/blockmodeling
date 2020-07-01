"IM" <-
function(res,which=1, drop=TRUE, ...){
  if(class(res)=="opt.more.par"){
    IM<-res$best[[which]]$IM
  } else IM<-res$IM
  if(drop)IM<-drop(IM)
  return(IM)
}

"EM" <-
function(res,which=1, drop=TRUE,...){
    if(class(res)=="opt.more.par"){
      EM<-res$best[[which]]$EM
    } else EM<-res$EM
    if(drop)EM<-drop(EM)
    return(EM)
}

