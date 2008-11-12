"savevector" <-
structure(function(v,filename){write(c(paste("*Vertices",length(v)), v), file = filename, ncolumns=1)}
, comment = "Save vector to file that can be read by Pajek")
