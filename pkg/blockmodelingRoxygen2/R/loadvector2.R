#' @rdname Pajek
#' 
#' @description \code{loadvector2} - The same as above, but adapted to be called within \code{loadpajek} - as a consequence not suited for reading clusters.
#' importFrom utils read.table


"loadvector2" <-
structure(function(filename){
  if(is.character(filename)) {file<-file(description=filename,open="r")
  }else file<-filename
  nn <-read.table(file=file,nrows=1,stringsAsFactors=FALSE)
  while(tolower(nn[1])!="*vertices") nn <-read.table(file=file,nrows=1,stringsAsFactors=FALSE)
  vv<-read.table(file=file,nrows=nn[[2]])
  if (dim(vv)[2]==1)
    vv<-vv[[1]]
  if(is.character(filename)) close(file)
  vv
}
, comment = "Load vector(s) from file that was produced by Pajek")
