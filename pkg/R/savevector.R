"savevector" <-
structure(function(v,filename){
if(grep(patt="w32",x=version["os"])){
	eol<-"\n"
}else{eol<-"\r\n"}
cat(paste(c(paste("*Vertices",length(v)), v),collapse=eol),file = filename)
}
, comment = "Save vector to file that can be read by Pajek")
