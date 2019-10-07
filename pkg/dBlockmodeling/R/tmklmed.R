tmklmed = function(A,RC,CC,TLIMIT) {
# This program runs two-mode KL-medians for an
# an RO x CO network matrix.
# INPUTS
#	A - an RO x CO two-mode network matrix
#     RC - the number of row clusters (1 < RC < RO)
#     CC - the number of row clusters (1 < CC < CO)
#     TLIMIT - a desired time limit
# OUTPUTS
#	objval - total number of inconsistencies
#	RP - an RO-dimensional vector of row cluser assignements
#	CP - an CO-dimensional vector of column cluser assignements
#     restarts - the number of restarts within the time limit
	RO = dim(A)[1]
	CO = dim(A)[2]
	GBEST = 0
      NREPS = 0
	GR <- matrix(0, nrow = RO, ncol = 1)
	GC <- matrix(0, nrow = CO, ncol = 1)
	res =.Fortran("tmklmed",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.integer(A),as.integer(GR),as.integer(GC),as.integer(GBEST),as.integer(NREPS))
	RP <- unlist(res[[7]])
	CP <- res[[8]]
	objval <- res[[9]]
	restarts <- res[[10]]
	return(list(RP=RP, CP=CP, objval=objval, restarts=restarts))
}
