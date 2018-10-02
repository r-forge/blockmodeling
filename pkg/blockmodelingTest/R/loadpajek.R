loadpajek<-function(filename){
	if(is.character(filename)) {file<-file(description=filename,open="r")
	}else file<-filename
	res<-list(Networks=list(),Partitions=list(),Vectors=list(),Permutation=list())
	nblanklines=0
	while(TRUE){
		line<-scan(file = file, nlines =1,what="char",quiet =TRUE, blank.lines.skip=FALSE)
		if(length(line)==0) {
			break
		}
        if (substr(line[1],start=1,stop=1)=="%") {
            print(paste(line,collapse=" "))
            next
        }
		if(line[1]=="") next
		if(sum(grep(pattern="^ *$",x=as.character(line))==1)) next
		if(tolower(tolower(line[1]))=="*matrix" || tolower(line[1])=="*network"){
			objName<-paste(line[-1],collapse=" ")
			if(tolower(line[1])=="*matrix"){
				readObj<-loadmatrix(file)
			}else readObj<-loadnetwork2(file, closeFile=FALSE)

			if(objName %in% names(res[["Networks"]])){
				i<-1
				while(TRUE){
					if(paste(objName,"Ver",i) %in% names(res[["Networks"]])) break
					i<-i+1
				}
				objName<-paste(objName,"Ver",i)
			}
			res[["Networks"]]<-c(res[["Networks"]],list(readObj))
			names(res[["Networks"]])[length(res[["Networks"]])]<-objName
		} else if(tolower(line[1])=="*vector" || tolower(line[1])=="*permutation" || tolower(line[1])=="*partition"){
			objName<-paste(line[-1],collapse=" ")
			readObj<-loadvector2(file)
			if(tolower(line[1])=="*vector"){
				if(objName %in% names(res[["Vectors"]])){
					i<-1
					while(TRUE){
						if(paste(objName,"Ver",i) %in% names(res[["Vectors"]])) break
						i<-i+1
					}
					objName<-paste(objName,"Ver",i)
				}
				res[["Vectors"]]<-c(res[["Vectors"]],list(readObj))
				names(res[["Vectors"]])[length(res[["Vectors"]])]<-objName
			} else if(tolower(line[1])=="*permutation"){
				if(objName %in% names(res[["Permutations"]])){
					i<-1
					while(TRUE){
						if(paste(objName,"Ver",i) %in% names(res[["Permutations"]])) break
						i<-i+1
					}
					objName<-paste(objName,"Ver",i)
				}
				res[["Permutations"]]<-c(res[["Permutations"]],list(readObj))
				names(res[["Permutations"]])[length(res[["Permutations"]])]<-objName
			} else if(tolower(line[1])=="*partition"){
				if(objName %in% names(res[["Partitions"]])){
					i<-1
					while(TRUE){
						if(paste(objName,"Ver",i) %in% names(res[["Partitions"]])) break
						i<-i+1
					}
					objName<-paste(objName,"Ver",i)
				}
				res[["Partitions"]]<-c(res[["Partitions"]],list(readObj))
				names(res[["Partitions"]])[length(res[["Partitions"]])]<-objName
			}
		}
	
	}
	return(res)
	close(file)
}
