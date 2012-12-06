.onAttach<-function(libname, pkgname){
    library(utils)
    cit<-citation(pkgname)
    txt<-paste(c(format(cit,"citation")),collapse="\n\n")
    packageStartupMessage(txt)
}