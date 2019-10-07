rhrsbt = function(A,C,TLIMIT) {
# This program runs nonnegative matrix factorization for
# an MxM asymmetric matrix. The main diagonal is ignored.
# INPUTS
#	A - an N x N signed network matrix
#     C - the number of cluster (1 < C < N)
#     TLIMIT - a desired time limit
# OUTPUTS
#	obj - the Doreian & Mrvar objective value
#	P - and N-dimensional vector of cluser assignements
#     restarts - the number of restarts within the time limit
	N = dim(A)[1]
	OBJVAL = 0
  NREPS = 0
	EBEST <- matrix(0, nrow = N, ncol = 1)
	res =.Fortran("rhrsbt",as.integer(N),as.integer(C),as.double(TLIMIT),as.double(OBJVAL),as.integer(A),as.integer(EBEST),as.integer(NREPS))
	P <- res[[6]]
	obj <- res[[4]]
	restarts <- res[[7]]
	return(list(P=P, obj=obj, restarts=restarts))
}
