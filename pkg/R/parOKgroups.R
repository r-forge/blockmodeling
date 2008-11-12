"parOKgroups" <-
function(clu,parOKaddParam = list(
k, 	#number of clusters by groups
groups)	#partition of units into exclusive groups
){
	isTRUE(all(cut(clu,c(0,cumsum(parOKaddParam$k)),labels =FALSE)==parOKaddParam$groups))
	
}
