"plot.mat" <- plotMat <- 
function(
    x=M, #x should be a matrix or similar object
    M=x, #M should be a matrix or similar object - both (x and M) are here to make the code compatible with generic plot and with older versions of plot.mat and possbily some other functions in the package
    clu=NULL,   #partition
    ylab="",
    xlab="",
    main=NULL,
    print.val=!length(table(M))<=2, #should the values be printed inside the cells
    print.0=FALSE,  #should the values equal to 0 be printed inside the cells, only used if 'print.val == TRUE'
    plot.legend=!print.val&&!length(table(M))<=2,   #should the legend for the colors be ploted
    print.legend.val="out", #where should the values for the legend be printed: 'out' - outside the cells (bellow), 'in' - inside the cells, 'both' - inside and outside the cells
    print.digits.legend=2,  #the number of digits that should appear in the legend
    print.digits.cells=2, #the number of digits that should appear in the cells (of the matrix and/or legend)
    print.cells.mf=NULL, #if not null, the above argument is igonred, the cell values are printed as the cell are multiplied by this factor and rounded
    outer.title=FALSE,  #should the title be printed on the 'inner' or 'outer' plot, default is 'inner' if legend is ploted and 'outer' otherwise
    title.line= ifelse(outer.title,-1.5,7), #the line (from the top) where the title should be printed
    mar= c(0.5, 7, 8.5, 0)+0.1, #A numerical vector of the form 'c(bottom, left, top, right)' which gives the lines of margin to be specified on the four sides of the plot. The default is 'c(5, 4, 4, 2) + 0.1'.
    cex.val="default",  #size of the values printed
    val.y.coor.cor = 0, #correction for centering the values in the sqares in y direction
    val.x.coor.cor = 0, #correction for centering the values in the sqares in x direction
    cex.legend=1,   #size of the text in the legend,
    legend.title="Legend",  #the title of the legend
    cex.axes="default", #size of the characters in axes, 'default' makes the cex so small that all categories can be printed
    print.axes.val=NULL,    #should the axes values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
    print.x.axis.val=!is.null(colnames(M)), #should the x axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
    print.y.axis.val=!is.null(rownames(M)), #should the y axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
    x.axis.val.pos = 1.01, #y coordiante of the x axis values
    y.axis.val.pos = -0.01,  #x coordiante of the y axis values
    cex.main=par()$cex.main,
    cex.lab=par()$cex.lab,
    yaxis.line=-1.5,    #the position of the y axis (the argument 'line')
    xaxis.line=-1,  #the position of the x axis (the argument 'line')
    legend.left=0.4,#how much left should the legend be from the matrix
    legend.up=0.03, #how much left should the legend be from the matrix
    legend.size=1/min(dim(M)),  #relative legend size
    legend.text.hor.pos=0.5,    #horizontal position of the legend text (bottom) - 0 = bottom, 0.5 = middle,...
    par.line.width = 3, #the width of the line that seperates the partitions
    par.line.col = "blue", #the color of the line that seperates the partitions
    IM.dens= NULL,
    IM= NULL,   #Image used for ploting (shaded lines)
    wnet=NULL,      #which net (if more) should be ploted - used if M is an array
    wIM=NULL,   #which IM (if more) should be used for ploting (defualt = wnet) - used if IM is an array
    use.IM=length(dim(IM))==length(dim(M))|!is.null(wIM),   #should IM be used for ploting?
    dens.leg=c(null=100, nul=100),
    blackdens=70,
    plotLines = FALSE, #Should the lines in the matrix be printed (best set to FALSE for larger networks)
    frameMatrix=TRUE, #Should the matrix be framed (if plotLines is FALSE)
    x0ParLine=-0.1, #x coordinates for lines between row clusters
    x1ParLine=1, #x coordinates for lines between row clusters
    y0ParLine=0, #y coordinates for lines between col clusters
    y1ParLine=1.1, #y coordinates for lines between col clusters
	colByUnits=NULL, #a vector (of 0s and 1s) indicating whether ties of a unit should be marked with a diferent (nonblack) color - only used for binary networks 
	colByRow=NULL, #a vector (of 0s and 1s) indicating whether outgoing ties of a unit should be marked with a different (nonblack) color - only used for binary networks
	colByCol=NULL, #a vector (of 0s and 1s) indicating whether incoming ties of a unit should be marked with a different (nonblack) color - only used for binary networks
    mulCol = 2,
    joinColOperator = "+",
    colTies=FALSE,
    maxValPlot=NULL, # maximal value used for determining the color of cells in the plot. This value and all higher (in absolute terms) will produce a pure black/red color
	printMultipliedMessage = TRUE, # shold mutiplication message be printed when values were the printed tie values are multiplied
	replaceNAdiagWith0=TRUE, #Should the diagonal with only NAs be replace by 0s?
	colLabels=FALSE, # Should the labels of units be colored. If FALSE, these are not collored, if TRUE, they are colored with colors of clusters as defined by palette. This can be aslo a vector of colors (or integers) for one-mode networks or a list of two such vectors for two-mode networks.
    ... #aditional arguments to plot.default
){
    old.mar<-par("mar")
    if(length(dim(IM))>length(dim(IM))&use.IM){
        if(is.null(wIM))wIM<-wnet
        if(is.null(wIM)) wIM<-1
        if(length(dim(IM))==3) {
          IM<-IM[wIM,,]
        } else{
          warning("IM will not be used for plotting. Cannot be sure how to extract the appropirate part!")
          use.IM<-FALSE
        }
    }
    tempClu<-clu
	

	
    if(length(dim(M))>2){
        if(!is.null(wnet)){
            relDim<-which.min(dim(M))
            if(relDim==1){
                M<-M[wnet,,]
            }else if(relDim==3){
                    M<-M[,,wnet]
            }else stop("More than 2 dimensions where relation dimension can not be determined")
        }else{
            plot.array(M = M,
                clu=tempClu,    #partition
                ylab=ylab,
                xlab=xlab,
                main.title=main,main.title.line=-2,
                print.val=print.val,    #should the values be printed inside the cells
                print.0=print.0,    #should the values equal to 0 be printed inside the cells, only used if 'print.val == TRUE'
                plot.legend=plot.legend,   #should the legend for the colors be ploted
                print.legend.val=print.legend.val,  #where should the values for the legend be printed: 'out' - outside the cells (bellow), 'in' - inside the cells, 'both' - inside and outside the cells
                print.digits.legend=print.digits.legend,    #the number of digits that should appear in the legend
                print.digits.cells=print.digits.cells, #the number of digits that should appear in the cells (of the matrix and/or legend)
                print.cells.mf=print.cells.mf, #if not null, the above argument is igonred, the cell values are printed as the cell are multiplied by this factor and rounded
                outer.title=outer.title,    #should the title be printed on the 'inner' or 'outer' plot, default is 'inner' if legend is ploted and 'outer' otherwise
                title.line= title.line, #the line (from the top) where the title should be printed
                mar= mar, #A numerical vector of the form 'c(bottom, left, top, right)' which gives the lines of margin to be specified on the four sides of the plot. The default is 'c(5, 4, 4, 2) + 0.1'.
                cex.val=cex.val,    #size of the values printed
                val.y.coor.cor = val.y.coor.cor, #correction for centering the values in the sqares in y direction
                val.x.coor.cor = val.x.coor.cor, #correction for centering the values in the sqares in x direction
                cex.legend=cex.legend,  #size of the text in the legend,
                legend.title=legend.title,  #the title of the legend
                cex.axes=cex.axes,  #size of the characters in axes, 'default' makes the cex so small that all categories can be printed
                print.axes.val=print.axes.val,  #should the axes values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
                print.x.axis.val=print.x.axis.val,  #should the x axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
                print.y.axis.val=print.y.axis.val,  #should the y axis values be printed, 'default' prints each axis if 'rownames' or 'colnames' is not 'NULL'
                x.axis.val.pos = x.axis.val.pos, #y coordiante of the x axis values
                y.axis.val.pos = y.axis.val.pos,  #x coordiante of the y axis values
                cex.main=cex.main,
                cex.lab=cex.lab,
                yaxis.line=yaxis.line,  #the position of the y axis (the argument 'line')
                xaxis.line=xaxis.line,  #the position of the x axis (the argument 'line')
                legend.left=legend.left,#how much left should the legend be from the matrix
                legend.up=legend.up,    #how much left should the legend be from the matrix
                legend.size=legend.size,    #relative legend size
                legend.text.hor.pos=legend.text.hor.pos,    #horizontal position of the legend text (bottom) - 0 = bottom, 0.5 = middle,...
                par.line.width = par.line.width , #the width of the line that seperates the partitions
                par.line.col = par.line.col, #the color of the line that seperates the partitions
                IM.dens= IM.dens,
                IM= IM, #Image used for ploting (shaded lines)
                wIM=wIM,    #which IM (if more) should be used for ploting (defualt = wnet) - used if IM is an array
                use.IM=use.IM,  #should IM be used for ploting?
                dens.leg=dens.leg,
                blackdens=blackdens,
                plotLines = plotLines,...
            )
            return(invisible(NULL))
        }
    }
    
    dm<-dim(M)

    if(class(M)!="matrix"&&class(M)!="mat"){
        pack<-attr(class(M),"package")
        if(!(is.null(pack))&&pack=="Matrix"){
            if(requireNamespace("Matrix")){
                M<-as.matrix(M)
            } else stop("The supplied object needs Matrix packege, but the package is not available (install it!!!).")
        } else {
            warning("Attempting to convert object of class ",class(M)," to class 'matrix'. Keep fingers crossed.")
            M<-as.matrix(M)
        }
    }
  
    if(replaceNAdiagWith0 & all(is.na(diag(M)))) diag(M)<-0

    if(is.null(main)){
        objName<-deparse(substitute(M))
        if(objName[1]=="x"){
                objName<-deparse(substitute(x))
        } 
        if(length(objName)>1)  objName=""
        main <- paste("Matrix",objName)
		if(nchar(main)>50) main<-substr(main,1,50)
    }
    #if(length(main)>26)

    if(is.logical(print.axes.val)){
        print.x.axis.val<-print.y.axis.val<-print.axes.val
    }


    #defining text on the axes if row or colnames do not exist
    if(is.null(rownames(M))){
        rownames(M)<-1:dm[1]
        }
    if(is.null(colnames(M))){
        colnames(M)<-1:dm[2]
    }

    if(!is.null(clu)){  #is any clustering provided, ordering of the matrix if 'TRUE'
      if(is.list(clu)){
        clu<-lapply(clu,function(x)as.integer(as.factor(x)))
        tmNclu<-sapply(clu,max)
        for(iMode in 2:length(tmNclu)){
          clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
        }
        unlistClu<-unlist(clu)
        if( all(length(unlistClu)==dm)) clu<-unlistClu
      }
        if(!is.list(clu)){
            tclu<-table(clu)
            or.c<-or.r<-order(clu)
            clu<-list(clu,clu)
            lines.col<-cumsum(tclu)[-length(tclu)]*1/dm[2]
            lines.row<-1-lines.col
        }else if(is.list(clu)&&length(clu)==2){
			if(!is.null(clu[[1]])){
				tclu.r<-table(clu[[1]])
				or.r<-order(clu[[1]])
				lines.row<- 1-cumsum(tclu.r)[-length(tclu.r)]*1/dm[1]
			} else{
				or.r<-1:dim(M)[1]
				lines.row<-NULL
			}
			if(!is.null(clu[[2]])){
				tclu.c<-table(clu[[2]])
				or.c<-order(clu[[2]])
				lines.col<-cumsum(tclu.c)[-length(tclu.c)]*1/dm[2]
			} else{
				or.c<-1:dim(M)[2]
				lines.col<-NULL
			}            
            
        } else stop("Networks with more that 2 modes (ways) must convert to 1-mode networks before it is sent to this function.")
        M<-M[or.r,or.c]
    clu<-lapply(clu,function(x)as.numeric(factor(x)))
    }



  if(is.null(IM.dens)){
    if(!is.null(IM)&use.IM){
      IM.dens<-matrix(-1,ncol=dim(IM)[2],nrow=dim(IM)[1])
      for(i in names(dens.leg)){
        IM.dens[IM==i]<- dens.leg[i]
      }
    }
  }


  if(!is.null(IM.dens)){
      dens<-matrix(-1,nrow=dm[1], ncol=dm[2])
    for(i in unique(clu[[1]])){
      for(j in unique(clu[[2]])){
        dens[clu[[1]]==i,clu[[2]]==j]<-IM.dens[i,j]
      }
    }
    dens<-dens[or.r,or.c]
  }

	if(length(cex.axes)==1) cex.axes<-c(cex.axes,cex.axes)
    if(cex.axes[1]=="default"){    #defining the size of text on the axes
        cex.y.axis<-min(15/dm[1],1)
    }else{
        cex.y.axis<-cex.axes[1]
    }
    if(cex.axes[2]=="default"){    #defining the size of text on the axes
        cex.x.axis<-min(15/dm[2],1)
    }else{
        cex.x.axis<-cex.axes[2]
    }

    #defining text on the axes
    yaxe<-rownames(M)
    xaxe<-colnames(M)


    ytop <- rep(x=(dm[1]:1)/dm[1],times=dm[2])  #definin the positions of rectangules
    ybottom<- ytop - 1/dm[1]
    xright <- rep(x=(1:dm[2])/dm[2],each=dm[1])
    xleft <- xright - 1/dm[2]
  
   if(all(M %in% c(0,1))){
#       browser()
            mulCol<-mulCol
        if(is.null(colByRow)&is.null(colByCol)) {
            colByRow<-colByCol<-colByUnits
        } else {
            if(is.null(colByRow)){
                colByRow<-rep(0, length(colByCol))
                mulCol<-1
            }
            if(is.null(colByCol)){
                colByCol<-rep(0, length(colByRow))
            }
            colByUnits<-TRUE
        }
        
        col<-M
        if(all(col %in% c(0,1))& (!is.null(colByUnits))){
            newCol<-outer(colByRow,colByCol*mulCol,FUN=joinColOperator)
            if(!is.null(clu)) newCol<-newCol[or.r,or.c]
            if(colTies){
                col[M>0]<-col[M>0]+newCol[M>0]
            }else{
                newCol[newCol>0]<-newCol[newCol>0]+1
                col[M==0]<-col[M==0]+newCol[M==0]
            }
        }        
    } else {
        aM<-abs(M)
        if(!is.null(maxValPlot)){
          aM[aM>maxValPlot]<-maxValPlot
        }
        max.aM<-max(aM)
        aMnorm<-as.vector(aM)/max.aM
        if(max.aM!=0){
            col<-grey(1-aMnorm)   #definin the color of rectangules
        }else col<-matrix(grey(1),nrow=dm[1],ncol=dm[2])
        col[M<0]<-paste("#FF",substr(col[M<0],start=4,stop=7),sep="")
    }

    asp<-dm[1]/dm[2]    #making sure that the cells are squares	

    par(mar=mar, xpd=NA)    #ploting
    plot.default(c(0,1),c(0,1),type="n",axes=FALSE,ann=FALSE,xaxs="i",asp=asp,...)
  if(is.null(IM.dens)||all(IM.dens==-1)){
    rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=col,cex.lab=cex.lab,border=if(plotLines)"black" else NA)
  }else{
    rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=col,cex.lab=cex.lab,density=dens,border=if(plotLines)"black" else NA)
  } 
  if(frameMatrix) rect(xleft=0, ybottom=0, xright=1, ytop=1, cex.lab=cex.lab,border="black")
    if(!is.null(clu)){  #ploting the lines between clusters
        if(!is.null(lines.row)) segments(x0=x0ParLine,x1=x1ParLine,y0=lines.row,y1=lines.row,col=par.line.col,lwd=par.line.width)
        if(!is.null(lines.col)) segments(y0=y0ParLine,y1=y1ParLine,x0=lines.col,x1=lines.col,col=par.line.col,lwd=par.line.width )
    }
	
	colYlabels <- colXlabels <- 1
	if((length(colLabels)==1)&&is.logical(colLabels)){
		if(colLabels){
			if(is.null(clu)){
				warning("clu not used!")
			} else {
				colYlabels <- clu[[1]]
				colXlabels <- clu[[2]]
			}
		} 
	} else{
		if(!is.list(colLabels))colLabels<-list(colLabels,colLabels)
		if(length(colLabels[[1]])==dm[1]){
			colYlabels<-colLabels[[1]]
		} else {
			warning("colLabels for first dimmension of wrong length, no colors will be used!")
		}
		if(length(colLabels[[2]])==dm[2]){
			colXlabels<-colLabels[[2]]
		} else {
			warning("colLabels for second dimmension of wrong length, no colors will be used!")
		}	
	}		
	if(!is.null(clu)){
		if(length(colXlabels)>1) colXlabels<-colXlabels[or.c]
		if(length(colYlabels)>1) colYlabels<-colYlabels[or.r]
	}
	
    if(print.y.axis.val) text(x=y.axis.val.pos, y = (dm[1]:1)/dm[1]-1/dm[1]/2 +val.y.coor.cor,labels = yaxe,cex=cex.y.axis,adj=1, col=colYlabels)
    if(print.x.axis.val) text(y=x.axis.val.pos, x = (1:dm[2])/dm[2]-1/dm[2]/2 +val.x.coor.cor, srt=90, labels = xaxe, cex=cex.x.axis,adj=0, , col=colXlabels)
    title(outer=outer.title,ylab=ylab,xlab=xlab,main=main, line=title.line,cex.main=cex.main)
    if(print.val){  #ploting the values in the cells if selected
        norm.val<-as.vector(M)/max(abs(M))
        aMnorm<-abs(norm.val)
        col.text<-1-round(aMnorm)

        if(!print.0) col.text[as.vector(M)==0]<-0

        if(length(table(col.text))==2) {
            col.labels<-c("white","black")
        } else col.labels<-c("white")

        col.text<-as.character(factor(col.text,labels=col.labels))
    if(!is.null(IM.dens)&&!all(IM.dens==-1)) col.text[col.text=="white"&dens>0&dens<blackdens]<-"black"

        col.text[col.text=="black"&norm.val<0]<-"red"
        if(!print.0) col.text[as.vector(M)==0]<-"transparent"

        maxM<-formatC(max(M),format="e")
        if(is.null(print.cells.mf)){
            if(all(trunc(M)==M)& max(M)<10^print.digits.cells){
                multi<-1
            }else{
                multi<-floor(log10(max(M)))
                multi<-(multi-(print.digits.cells - 1))*(-1)
                multi<-10^multi
            }
        }else multi <- print.cells.mf
        M.plot<-round(M*multi)

        text(x=(xleft+xright)/2+val.x.coor.cor,y=(ytop+ybottom)/2+val.y.coor.cor, labels=as.vector(M.plot),col=col.text,cex=ifelse(cex.val=="default",min(10/max(dm),1),cex.val))
        if(multi!=1 & printMultipliedMessage) mtext(text=paste("* all values in cells were multiplied by ",multi,sep=""),side=1, line=-0.7,cex=0.70)
    }

    if(plot.legend){    #ploting the legend if selected
        if(asp>=1){
            xright.legend<- -legend.left
            xleft.legend <- xright.legend - 1*legend.size*asp
            ybottom.legend <- 1+(4:0)*legend.size+ legend.up
            ytop.legend <- ybottom.legend + 1*legend.size
        }else{
            xright.legend<- -legend.left
            xleft.legend <- xright.legend - 1*legend.size
            ybottom.legend <- 1+(4:0)*legend.size*asp+ legend.up
            ytop.legend <- ybottom.legend + 1*legend.size*asp
        }
        col.legend<-gray(4:0/4)
        rect(xleft=xleft.legend, ybottom=ybottom.legend, xright=xright.legend, ytop=ytop.legend, col=col.legend)
        if(print.legend.val=="out"|print.legend.val=="both") text(x=xright.legend + 1/20,y= (ytop.legend+ybottom.legend)/2, labels=formatC(0:4/4*max(M), digits = print.digits.legend,format="g"),adj=0,cex=cex.legend)
        text(x=xleft.legend,y=ytop.legend[1] + legend.size/asp/2+0.02, labels=legend.title,font=2,cex=cex.legend,adj=0)

        if(print.legend.val=="in"|print.legend.val=="both"){
            col.text.legend<-round(4:0/4)
            if(!print.0) col.text.legend[1]<-0
            col.text.legend<-as.character(factor(col.text.legend,labels=c("white","black")))
            if(!print.val){
                if(is.null(print.cells.mf)){
                    if(all(trunc(M)==M)& max(M)<10^print.digits.cells){
                        multi<-1
                    }else{
                        multi<-floor(log10(max(M)))
                        multi<-(multi-(print.digits.cells - 1))*(-1)
                        multi<-10^multi
                    }
                }else multi <- print.cells.mf
                maxM<-round(max(M)*multi)
            } else maxM<-max(M.plot)
            text(x=(xleft.legend+xright.legend)/2,y=(ytop.legend+ybottom.legend)/2, labels=round(0:4/4*maxM),col=col.text.legend,cex=cex.legend)
        }
    }

    par(mar=old.mar)
}


"plot.array" <- plotArray <- 
function(
    x=M, #x should be a matrix or similar object
    M=x, #M should be a matrix or similar object - both (x and M) are here to make the code compatible with generic plot and with older versions of plot.mat and possbily some other functions in the package
    IM=NULL, #the image to be used for plotting
    ...,    #aditional arguments to plot.mat
    main.title=NULL,main.title.line=-2,mfrow=NULL
){
    if(is.null(main.title)){
        objName<-deparse(substitute(M))
        if(objName=="x")objName<-deparse(substitute(x))
        main.title <- paste("Matrix",objName)
		if(nchar(main.title)>50) main.title<-substr(main.title,1,50)
    }
    dM<-dim(M)
    relDim<-which.min(dM)
    nDim<-dM[relDim]
    if(is.null(mfrow)|(prod(mfrow)<nDim)){
        if(nDim<4){
            mfrow<-c(1,nDim)
        } else if(nDim<6){
            mfrow<-c(2,ceiling(nDim/2))
        } else{
            nr<-round(sqrt(nDim/6)*2); nc<-ceiling(nDim/nr)
            mfrow<-c(nr,nc)
        }           
    }
    par.def<-par(no.readonly = TRUE) 
    par(mfrow=mfrow)
    
    relNames<-dimnames(M)[[relDim]]
    if(is.null(relNames)) relNames<-1:nDim
    for(i in 1:nDim){
    #for(iName in relNames) 
        iName<-relNames[i]
        if(relDim==1){
            plot.mat(M[iName,,],main=iName, IM=IM[i,,],...)
        } else if(relDim==3) plot.mat(M[,,iName],main=iName, IM=IM[i,,],...)
    }
    
    title(main=main.title,outer=TRUE,line=main.title.line)
    par(par.def)
}


