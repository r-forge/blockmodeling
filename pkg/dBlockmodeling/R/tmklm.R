tmklm = function(A,RC,CC,TLIMIT) {
# This program runs two-mode K-means for an
# an RO x CO network matrix.
# INPUTS
#	A - an RO x CO two-mode network matrix
#     RC - the number of row clusters (1 < RC < RO)
#     CC - the number of row clusters (1 < CC < CO)
#     TLIMIT - a desired time limit
# OUTPUTS
#	vaf - the variance accounted
#	RP - an RO-dimensional vector of row cluser assignements
#	CP - an CO-dimensional vector of column cluser assignements
#     restarts - the number of restarts within the time limit
	RO = dim(A)[1]
	CO = dim(A)[2]
	VAF = 0
  NREPS = 0
	RBEST <- matrix(0, nrow = RO, ncol = 1)
	CBEST <- matrix(0, nrow = CO, ncol = 1)
	res =.Fortran("tmklm",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.double(A),as.integer(RBEST),as.integer(CBEST),as.double(VAF),as.integer(NREPS))
#	res =.Fortran("tmklm",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.double(A))
	RP <- res[[7]]
	CP <- res[[8]]
	vaf <- res[[9]]
	restarts <- res[[10]]
	return(list(RP=RP, CP=CP, vaf=vaf, restarts=restarts))
}
