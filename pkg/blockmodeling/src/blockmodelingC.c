/*

WARNINGS: rdo and cdo blocks added only to binary and valued blockmodeling - these blocks return Inf for homogeneity blockmodeling!



This is an implementation of
Generalized blockmodeling of valued (and binary) networks
in C to be called from R.

The structure is as follows:
- main functions that are linked to R. These must include all functions for optimizing partitions
- a function that computes (or updates) the criterion function of a partition and a blockmodel
- functions for computing errors/inconsistencies for individual block types for all types of blockmodeling (hom,val,bin)
- functions for computing row/column summaries for regular-like blocks
- functions for computing measure of variation for homogeneity blockmodeling


The implementation must support (when final):
- several types of generalized blockmodeling, with the possibility to extend it (easily) with new types of blockmodeling or/and types of blocks
- multirelational networks (with the possibility to have different or same images for all networks
- at least one and two mode networks (preferably efficiently) and also 3-mode networks
- efficient computations for symmetric networks
- symmetrical block type (not yet implemented)
- pre-specified blockmodeling - for valued blockmodeling this also means possibility to pre-specify the value from the which the deviations are computed (hom) or the value of parameter m by blocks (val)

If possible, also:
- possibility of different methods for "searching" for the partition - ways of optimizing - eg. not only local search, but also genetic algorithm, tabu search, etc.
- possibility to specify what kind of partitions are allowed (minimal/maximal group size, etc.)


TODO:
- allow penalties by relations (already implemented), by block types and by "positions". This could be implemented in C by just one 4d weighting array that could be in R computed (if desired) from those separate weighting schemes.
*/


#include <stdio.h>
#include <R.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Change these when you add new functions */
#define nRegFun 3
#define nHomFun 3
#define nBlockTypes 11
#define nApproaches 3
/* #define MaxNumOfDiffBlockTypes 10 */

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

const double eps = 0.001;

double ss(double *px, int n, double preSpecVal);
double ssP(double *px, int n, double preSpecVal);
double ss0(double *px, int n, double preSpecVal);
double ssPmin(double *px, int n, double preSpecVal);

double ad(double *px, int n, double preSpecVal);
double adP(double *px, int n, double preSpecVal);
double ad0(double *px, int n, double preSpecVal);

int randomInt(int n);

/* A function with returns a random number on the interval [0, n-1]*/
int randomInt(int n)
{
	int r;
	r = (int) (unif_rand()*n);
	return(r);
}



void randomCall(int *n, int *r);

void randomCall(int *n, int *r){
	GetRNGstate();	/* Get .Random.seed from R */
	*r = randomInt(*n);
	PutRNGstate(); /* Write .Random.seed in R */
}

/* Definition of an array of pointers to a function for computing some measure of variance*/
double (*phom[nHomFun][4])(double *px, int n, double preSpecVal);


/* A function for binary log-likelihood with probability equal to density*/
double bll(double *px, int n, double preSpecVal)
{
	double p=0;
	double res=0;
	int i;

	for(i=0;i<n;i++){
		p += px[i];
	}
	p = p/n;
	if(p < eps) {p=eps;}
	if(p > (1- eps)) {p = 1-eps;}
	
	
	for(i=0;i<n;i++){
		res += px[i]*log(p)+ (1-px[i])*log(1-p);
	}
	return(-res);
}


/* A function for binary log-likelihood with probability equal to preSpecVal */
double bllP(double *px, int n, double preSpecVal)
{
	double p=preSpecVal;
	double res=0;
	int i;

	if(p < eps) {p=eps;}
	if(p > (1- eps)) {p = 1-eps;}
	
	
	for(i=0;i<n;i++){
		res += px[i]*log(p)+ (1-px[i])*log(1-p);
	}
	return(-res);
}


/* A function for binary log-likelihood with probability equal to density*/
double bllPmin(double *px, int n, double preSpecVal)
{
	double p=0;
	double res=0;
	int i;

	for(i=0;i<n;i++){
		p += px[i];
	}
	p = p/n;
	
	if(preSpecVal > p){
		p=preSpecVal;
	}
	
	if(p < eps) {p=eps;}
	if(p > (1- eps)) {p = 1-eps;}
	
	
	for(i=0;i<n;i++){
		res += px[i]*log(p)+ (1-px[i])*log(1-p);
	}
	return(-res);
}


/* A function for binary log-likelihood with probability equal to 0+eps */
double bll0(double *px, int n, double preSpecVal)
{
	double p=eps;
	double res=0;
	int i;
	
	for(i=0;i<n;i++){
		res += px[i]*log(p)+ (1-px[i])*log(1-p);
	}
	return(-res);
}


/* A function for computing sum of squares deviations from the mean*/
double ss(double *px, int n, double preSpecVal)
{
	double ssx=0;
	double sumx=0;
	int i;

	for(i=0;i<n;i++){
		sumx += px[i];
		ssx += px[i]*px[i];
	}
	return(ssx - sumx*sumx/n);
}

/* A function for computing sum of squares deviations from preSpecVal given value*/
double ssP(double *px, int n, double preSpecVal)
{
	double ssx=0;
	double tmp;
	int i;
	for(i=0;i<n;i++){
		tmp= px[i] - preSpecVal;
		ssx += tmp * tmp;
	}
	return(ssx);
}

/* A function for computing sum of squares deviations from preSpecVal given value or higher*/
double ssPmin(double *px, int n, double preSpecVal)
{
	double ssxP=0;
	double sumx=0;
	double tmp;
	double m;
	int i;
	for(i=0;i<n;i++){
		sumx += px[i];
	}
	m=sumx/n;
	if(preSpecVal > m){
		m=preSpecVal;
	}
	for(i=0;i<n;i++){
		tmp= px[i] - m;
		ssxP += tmp * tmp;
	}
	return(ssxP);
}

/* A function for computing sum of squares deviations from 0 */
double ss0(double *px, int n, double preSpecVal)
{
	double ssx=0;
	int i;
	for(i=0;i<n;i++){
		ssx += px[i] * px[i];
	}
	return(ssx);
}

/* function for compassion of two doubles "double" - copied form GNU C Library manual */
int cmp (const void *a, const void *b) {
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

	
/* A function for computing sum of absolute deviations from the median*/
double ad(double *px, int n, double preSpecVal){
	/*int cmp(double *x, double *y){
		if(*x>*y) return(1);
		if(*y>*x) return((-1));
		return(0);
	}*/




	double med, sad = 0;
	int i;
	qsort(px,n,sizeof(double), cmp);
	if((n%2)==0){
		med = ( px[n/2-1] + px[n/2])/2.0;
	}else{
		med = px[(n)/2];
	}
	for(i=0; i<(n/2) ;i++){
		sad += med - px[i];
	}

	for(i=n/2;i<(n);i++){
		sad += px[i] - med;
	}
	return(sad);
}


/* A function for computing sum of absolute deviations from a given value*/
double adP(double *px, int n, double preSpecVal){

	double sad = 0;
	int i;
	for(i=0; i<(n) ;i++){
		sad += px[i]>preSpecVal ? (px[i]- preSpecVal): (preSpecVal - px[i]);
	}
	return(sad);
}


/* A function for computing sum of absolute deviations from a given value*/
double ad0(double *px, int n, double preSpecVal){

	double sad = 0;
	int i;
	for(i=0; i<(n) ;i++){
		sad += px[i]>0 ? (px[i]): ( - px[i]);
	}
	return(sad);
}


double adPmin(double *px, int n, double preSpecVal){
	double med, sad = 0;
	int i;
	qsort(px,n,sizeof(double), cmp);
	if((n%2)==0){
		med = ( px[n/2-1] + px[n/2])/2.0;
	}else{
		med = px[(n)/2];
	}
	if(preSpecVal > med){
		med=preSpecVal;
	}

	for(i=0; i<(n/2) ;i++){
		sad += med - px[i];
	}

	for(i=n/2;i<(n);i++){
		sad += px[i] - med;
	}
	return(sad);
}


/* Definition of an array of pointers to a function for computing some summary measure*/
double (*pregFuns[nRegFun])(double *px, int n);

double maxv(double *px, int n){
	double res=-INFINITY;

	for(int i = 0;i<n;i++){
		res=max(res,px[i]);
	}
	return(res);
}

double sumv(double *px, int n){
	double res=0;

	for(int i = 0;i<n;i++){
		res +=px[i];
	}
	return(res);
}

double meanv(double *px, int n){
	double res=0;

	for(int i = 0;i<n;i++){
		res +=px[i];
	}
	return(res/n);
}



/* a function for computing error of  the dnc (do not care) block for all blockmodeling types -> always returns 0*/
double doNotCare(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	return(0.0);
}

/* a function for computing error of the block that always returns Inf*/
double infBlock(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	return(INFINITY);
}


/* a function for computing error of  the regular block - binary blockmodeling*/
double binReg(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of rows in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of rows in the block
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double *prs;
	double *pcs;
	prs = (double *) malloc(nrb*sizeof(double));
	pcs = (double *) malloc(ncb*sizeof(double));
	for(int i = 0; i<nrb;i++){
		prs[i]=0;
	}
	int nnr=0;
	int nnc=0;

	for(int j = 0; j<ncb; j++){
		pcs[j]=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			prs[i] += pM[ind2d+pRowInd[i]];
			pcs[j] += pM[ind2d+pRowInd[i]];
		}
		nnc += (pcs[j]>0);
	}
	for(int i = 0; i<nrb;i++){
		nnr += (prs[i]>0);
	}

	free(prs);
	free(pcs);
	if(*pmulReg==1) {
		return((nrb-nnr)*ncb + (ncb-nnc)*nnr);
	} else {
		return((nrb-nnr) + (ncb-nnc));
	}
}


/* a function for computing error of  the column-regular block - binary blockmodeling*/
double binCre(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	int nnc=0;

	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcs += pM[ind2d+pRowInd[i]];
		}
		nnc += (pcs>0);
	}
	
	if(*pmulReg==1) {
		return((ncb-nnc)*nrb);
	} else {
		return(ncb-nnc);
	}
}


/* a function for computing error of  the row-regular block - binary blockmodeling*/
double binRre(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	double prs=0;
	int nnr=0;

	for(int i = 0; i<nrb; i++){
		prs=0;
		for(int j = 0; j<ncb; j++){
			prs += pM[baseInd+nc*pColInd[j]+pRowInd[i]];
		}
		nnr += (prs>0);
	}
	
	if(*pmulReg==1) {
		return((nrb-nnr)*ncb);
	} else {
		return((nrb-nnr));
	}
}



/* a function for computing error of  the row-dominant block - binary blockmodeling*/
double binRdo(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	double prs=0;
	double mrs=0;

	for(int i = 0; i<nrb; i++){
		prs=0;
		for(int j = 0; j<ncb; j++){
			prs += pM[baseInd+nc*pColInd[j]+pRowInd[i]];
		}
		if(mrs<prs){mrs = prs; }
	}
	
	if(*pmulReg==1) {
		return((nrb-mrs)*ncb);
	} else {
		return((nrb-mrs));
	}
}


/* a function for computing error of  the row-dominant - binary blockmodeling - diagonal*/
double binRdoDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	int baseInd=relN*nr*nc;
	double prs=0;
	double mrs=0;
	double diagRes=0;

	for(int i = 0; i<nrb; i++){
		prs=0;
		diagRes += pM[baseInd+ nc*pColInd[i]+pRowInd[i]];
		for(int j = 0; j<ncb; j++){
			prs += pM[baseInd+nc*pColInd[j]+pRowInd[i]];
		}
		if(mrs<prs){mrs = prs; }
	}
	
	double res = nrb-mrs;
	if(diagRes==0) {res = res - 1;}
	if(*pmulReg==1) {
		return((res)*ncb);
	} else {
		return((res));
	}
}

/* a function for computing error of  the row-dominant block - binary blockmodeling - ignore diagonal*/
double binRdoIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	double prs=0;
	double mrs=0;

	for(int i = 0; i<nrb; i++){
		prs=0;
		for(int j = 0; j<ncb; j++){
			if(i!=j){
				prs += pM[baseInd+nc*pColInd[j]+pRowInd[i]];
			}
		}
		if(mrs<prs){mrs = prs; }
	}
	
	if(*pmulReg==1) {
		return((nrb-mrs-1)*ncb);
	} else {
		return((nrb-mrs-1));
	}
}


/* a function for computing error of  the row-dominant block - binary blockmodeling*/
double binCdo(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	double mcs=0;

	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcs += pM[ind2d+pRowInd[i]];
		}
		if(mcs<pcs){mcs = pcs; }
	}
	
	if(*pmulReg==1) {
		return((ncb-mcs)*nrb);
	} else {
		return((ncb-mcs));
	}
}


/* a function for computing error of  the row-dominant - binary blockmodeling - diagonal*/
double binCdoDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	double mcs=0;
	double diagRes=0;
	
	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		diagRes += pM[ind2d+pRowInd[j]];
		for(int i = 0; i<nrb; i++){
			pcs += pM[ind2d+pRowInd[i]];
		}
		if(mcs<pcs){mcs = pcs; }
	}

	double res = ncb-mcs;
	if (diagRes==0) {res = res - 1;}
	
	if(*pmulReg==1) {
		return((res)*nrb);
	} else {
		return((res));
	}
	
}

/* a function for computing error of  the row-dominant block - binary blockmodeling - ignore diagonal*/
double binCdoIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	double mcs=0;

	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			if(i!=j){
				pcs += pM[ind2d+pRowInd[i]];
			}
		}
		if(mcs<pcs){mcs = pcs; }
	}
	
	if(*pmulReg==1) {
		return((ncb-mcs-1)*nrb);
	} else {
		return((ncb-mcs-1));
	}
}


/* a function for computing error of  the row-functional block - binary blockmodeling*/
double binRfn(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	int baseInd=relN*nr*nc;
	double prs = 0;
	int nnr=0;
	double st=0;

	for(int i = 0; i<nrb; i++){
		prs=0;
		for(int j = 0; j<ncb; j++){
			prs += pM[baseInd+nc*pColInd[j]+pRowInd[i]];
		}
		nnr += (prs>0);
		st += prs;
	}

	
	if(*pmulReg==1) {
		return(st - nnr + (nrb-nnr)*ncb);
	} else {
		return(st - nnr + (nrb-nnr));
	}
}

/* a function for computing error of  the column-functional block - binary blockmodeling*/
double binCfn(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs = 0;
	double st = 0;
	int nnc=0;

	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcs += pM[ind2d+pRowInd[i]];
		}
		nnc += (pcs>0);
		st += pcs;
	}
	
	if(*pmulReg==1) {
		return(st - nnc + (ncb-nnc)*nrb);
	} else {
		return(st - nnc + (ncb-nnc));
	}
}


/* a function for computing error of  the complete block - binary blockmodeling*/
double binCom(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			res += pM[ind2d+pRowInd[i]];
		}
	}
	return(ncb*nrb - res);
}

/* a function for computing error of  the complete block - binary blockmodeling - diagonal*/
double binComDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	double diagRes=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	if(nrb==1){
		return(0.0);
	} else {
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			diagRes += pM[ind2d+pRowInd[j]];
			for(int i = (j + 1); i<nrb; i++){
				res += pM[ind2d+pRowInd[i]];
				res += pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
			}
		}
		return(ncb*(nrb - 1) - res + min(diagRes,(nrb-diagRes)));
	}
}



/* a function for computing error of  the complete block - binary blockmodeling - ignore diagonal*/
double binComIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	if(nrb==1){
		return(0);
	} else {
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			for(int i = (j + 1); i<nrb; i++){
				res += pM[ind2d+pRowInd[i]];
				res += pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
			}
		}
		return(ncb*(nrb - 1) - res );
	}
}




/* a function for computing error of  the null block - binary blockmodeling*/


double binNul(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0.0;
	int baseInd=relN*nr*nc;
	int ind2d;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			res += pM[ind2d+pRowInd[i]];
/*Rprintf("i = %i, j = %i, relNum = %i, tie value= %.2f\n", i,j,relN,pM[ind2d+pRowInd[i]]);*/
/*Rprintf("temp sum = %.2f\n", res);*/
		}
	}

/*Rprintf("binBlockErr = %.2f\n", res);*/
	return(res);
}

/* a function for computing error of  the null block - binary blockmodeling - diagonal*/
double binNulDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	double diagRes=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	if(nrb==1){
		return(0.0);
	} else {
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			diagRes += pM[ind2d+pRowInd[j]];
			for(int i = (j + 1); i<nrb; i++){
				res += pM[ind2d+pRowInd[i]];
				res += pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
			}
		}
		return(res + min(diagRes,(nrb-diagRes)));
	}
}



/* a function for computing error of  the null block - binary blockmodeling - diagonal*/
double binNulIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	if(nrb==1){
		return(0);
	} else {
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			for(int i = (j + 1); i<nrb; i++){
				res += pM[ind2d+pRowInd[i]];
				res += pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
			}
		}
		return(res);
	}
}





/* a function for computing error of  the regular block - homogeneity blockmodeling*/
double homReg(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	double rowErr;
	double colErr;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *prowArr;
	double *pcolArr;
	prowArr = (double *) malloc(nrb*ncb*sizeof(double));
	pcolArr = (double *) malloc(ncb*nrb*sizeof(double));

	double *prowStats;
	double *pcolStats;
	prowStats = (double *) malloc(nrb*sizeof(double));
	pcolStats = (double *) malloc(ncb*sizeof(double));



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcolArr[j*nrb + i] = pM[ind2d+pRowInd[i]];
			prowArr[i*ncb + j] = pM[ind2d+pRowInd[i]];
		}
		pcolStats[j]=pregFuns[regFun]((pcolArr + j*nrb),nrb);
	}
	for(int i = 0; i<nrb;i++){
		prowStats[i]=pregFuns[regFun]((prowArr + i*ncb),ncb);
	}
	free(prowArr);
	free(pcolArr);

	rowErr=phom[homFun][usePreSpecVal](prowStats,nrb,preSpecVal);
	colErr=phom[homFun][usePreSpecVal](pcolStats,ncb,preSpecVal);
	free(prowStats);
	free(pcolStats);

	
	if(*pmulReg==1) {
		return(max(rowErr*ncb,colErr*nrb));
	} else {
		return(max(rowErr,colErr));
	}
}


/* a function for computing error of  the column-regular block - homogeneity blockmodeling*/
double homCre(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	double colErr;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *pcolArr;
	pcolArr = (double *) malloc(ncb*nrb*sizeof(double));

	double *pcolStats;
	pcolStats = (double *) malloc(ncb*sizeof(double));



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcolArr[j*nrb + i] = pM[ind2d+pRowInd[i]];
		}
		pcolStats[j]=pregFuns[regFun]((pcolArr + j*nrb),nrb);
	}
	free(pcolArr);

	colErr=phom[homFun][usePreSpecVal](pcolStats,ncb,preSpecVal);
	free(pcolStats);

	
	if(*pmulReg==1) {
		return(colErr*nrb);
	} else {
		return(colErr);
	}
}



/* a function for computing error of  the row-regular block - homogeneity blockmodeling*/
double homRre(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	double rowErr;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *prowArr;
	prowArr = (double *) malloc(nrb*ncb*sizeof(double));

	double *prowStats;
	prowStats = (double *) malloc(nrb*sizeof(double));



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			prowArr[i*ncb + j] = pM[ind2d+pRowInd[i]];
		}
	}
	for(int i = 0; i<nrb;i++){
		prowStats[i]=pregFuns[regFun]((prowArr + i*ncb),ncb);
	}
	free(prowArr);

	rowErr=phom[homFun][usePreSpecVal](prowStats,nrb,preSpecVal);
	free(prowStats);

	
	if(*pmulReg==1) {
		return(rowErr*ncb);
	} else {
		return(rowErr);
	}
}


/* a function for computing error of  the row-functional block - homogeneity blockmodeling - experimental*/
double homRfn(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	double *px;
	double rowErr;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *prowArr;
	prowArr = (double *) malloc(nrb*ncb*sizeof(double));

	double *prowStats;
	prowStats = (double *) malloc(nrb*sizeof(double));
	
	int k=0;
	px = (double *) malloc(nrb*ncb*sizeof(double));
	
	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			prowArr[i*ncb + j] = pM[ind2d+pRowInd[i]];
			px[k]=pM[ind2d+pRowInd[i]];
			k++;
		}
	}
	for(int i = 0; i<nrb;i++){
		prowStats[i]=pregFuns[0]((prowArr + i*ncb),ncb);
		/* 0 in pregFuns[0] stands for maximum */
	}
	free(prowArr);

	rowErr=phom[homFun][usePreSpecVal](prowStats,nrb,preSpecVal);
	double nulErr=phom[homFun][2](px,nrb*ncb,0) - phom[homFun][2](prowStats,nrb,0);
	free(prowStats);	
	free(px);
	
	
	if(*pmulReg==1) {
		return(rowErr*ncb + nulErr);
	} else {
		return(rowErr + nulErr);
	}
}



/* a function for computing error of  the column-functional block - homogeneity blockmodeling - experimental*/
double homCfn(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	double *px;
	double colErr;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *pcolArr;
	pcolArr = (double *) malloc(ncb*nrb*sizeof(double));

	double *pcolStats;
	pcolStats = (double *) malloc(ncb*sizeof(double));

	int k=0;
	px = (double *) malloc(nrb*ncb*sizeof(double));



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcolArr[j*nrb + i] = pM[ind2d+pRowInd[i]];
			px[k]=pM[ind2d+pRowInd[i]];;
			k++;
		}
		pcolStats[j]=pregFuns[0]((pcolArr + j*nrb),nrb);
		/* 0 in pregFuns[0] stands for maximum */
	}
	free(pcolArr);

	colErr=phom[homFun][usePreSpecVal](pcolStats,ncb,preSpecVal);
	double nulErr=phom[homFun][2](px,nrb*ncb,0) - phom[homFun][2](pcolStats,ncb,0);
	free(pcolStats);
	free(px);

	
	if(*pmulReg==1) {
		return(colErr*nrb + nulErr);
	} else {
		return(colErr + nulErr);
	}
}





/* a function for computing error of  the complete block - homogeneity blockmodeling*/
double homCom(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

/*Rprintf("homCom - start \n");*/
	double *px;
	int baseInd=relN*nr*nc;
	int ind2d;
	int k=0;
	px = (double *) malloc(nrb*ncb*sizeof(double));
	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			px[k]=pM[ind2d+pRowInd[i]];
			k++;
		}
	}
/*Rprintf("homCom - end - 1  \n");	*/
	double res=phom[homFun][usePreSpecVal](px,nrb*ncb,preSpecVal);
	free(px);
/*Rprintf("homCom - end \n");*/
	return(res);
}



/* a function for computing error of  the complete block - homogeneity blockmodeling - diagonal*/
double homComDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

/*Rprintf("homComDiag - start \n");*/
	double *px;
	double *pdiag;
	int baseInd=relN*nr*nc;
	int ind2d;
	int k=0;

	if(nrb==1){
		return(0.0);
	} else {
		px = (double *) malloc(nrb*(ncb-1)*sizeof(double));
		pdiag = (double *) malloc(nrb*sizeof(double));
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
	/*Rprintf("indexDiag = %i \n", ind2d+pRowInd[j]);*/
			pdiag[j]=pM[ind2d+pRowInd[j]];
			for(int i = (j + 1); i<nrb; i++){
				px[k]=pM[ind2d+pRowInd[i]];
				k++;
				px[k]=pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
				k++;
			}
		}

	/*Rprintf("homComDiag - 1 \n");*/
	/*Rprintf("px: ");*/
	/*for( int i=0;i<nrb*(ncb-1);i++){*/
	/*	Rprintf("%.4f ", px[i]);*/
	/*}*/
	/*Rprintf("\n");*/

	/*Rprintf("pdiag: ");*/
	/*for( int i=0;i<nrb;i++){*/
	/*	Rprintf("%.4f ", pdiag[i]);*/
	/*}*/
	/*Rprintf("\n");*/
	/**/
	/*Rprintf("homFun = %i\n", homFun);*/
	/*Rprintf("usePreSpecVal = %i\n", usePreSpecVal);*/
	/*Rprintf("preSpecVal = %.4f\n", preSpecVal);*/

		double res=phom[homFun][usePreSpecVal](px,nrb*(ncb-1),preSpecVal)+phom[homFun][0](pdiag,nrb,0);
	/*Rprintf("homComDiag - 3 \n");	*/
		free(px);
	/*Rprintf("homComDiag - 4 \n");*/
		free(pdiag);
	/*Rprintf("homCOmDiag - end \n");*/
		return(res);
	}
}


/* a function for computing error of  the complete block - homogeneity blockmodeling - diagonal*/
double homComIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double *px;
	int baseInd=relN*nr*nc;
	int ind2d;
	int k=0;

	if(nrb==1){
		return(0.0);
	} else {
		px = (double *) malloc(nrb*(ncb-1)*sizeof(double));
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			for(int i = (j + 1); i<nrb; i++){
				px[k]=pM[ind2d+pRowInd[i]];
				k++;
				px[k]=pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
				k++;
			}
		}

		double res=phom[homFun][usePreSpecVal](px,nrb*(ncb-1),preSpecVal);
		free(px);
		return(res);
	}
}



/* a function for computing error of  the complete block - homogeneity blockmodeling*/
double homNul(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double *px;
	int baseInd=relN*nr*nc;
	int ind2d;
	int k=0;
	px = (double *) malloc(nrb*ncb*sizeof(double));
	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			px[k]=pM[ind2d+pRowInd[i]];
			k++;
		}
	}
	double res=phom[homFun][2](px,nrb*ncb,0);
	free(px);
	return(res);
}




/* a function for computing error of  the complete block - homogeneity blockmodeling - diagonal*/
double homNulDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double *px;
	double *pdiag;
	int baseInd=relN*nr*nc;
	int ind2d;
	int k=0;

	if(nrb==1){
		return(0.0);
	} else {
		px = (double *) malloc(nrb*(ncb-1)*sizeof(double));
		pdiag = (double *) malloc(nrb*sizeof(double));
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			pdiag[j]=pM[ind2d+pRowInd[j]];
			for(int i = (j + 1); i<nrb; i++){
				px[k]=pM[ind2d+pRowInd[i]];
				k++;
				px[k]=pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
				k++;
			}
		}
		double res=phom[homFun][2](px,nrb*(ncb-1),0)+phom[homFun][0](pdiag,nrb,0);
		free(px);
		free(pdiag);
		return(res);
	}
}


/* a function for computing error of  the complete block - homogeneity blockmodeling - diagonal*/
double homNulIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double *px;
	int baseInd=relN*nr*nc;
	int ind2d;
	int k=0;

	if(nrb==1){
		return(0.0);
	} else {
		px = (double *) malloc(nrb*(ncb-1)*sizeof(double));
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			for(int i = (j + 1); i<nrb; i++){
				px[k]=pM[ind2d+pRowInd[i]];
				k++;
				px[k]=pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
				k++;
			}
		}

		double res=phom[homFun][2](px,nrb*(ncb-1),0);
		free(px);
		return(res);
	}
}



/* a function for computing error of  the regular block - valued blockmodeling*/
double valReg(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *prowArr;
	double *pcolArr;
	prowArr = (double *) malloc(nrb*ncb*sizeof(double));
	pcolArr = (double *) malloc(ncb*nrb*sizeof(double));

	double *prowStats;
	double *pcolStats;
	prowStats = (double *) malloc(nrb*sizeof(double));
	pcolStats = (double *) malloc(ncb*sizeof(double));



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcolArr[j*nrb + i] = pM[ind2d+pRowInd[i]];
			prowArr[i*ncb + j] = pM[ind2d+pRowInd[i]];
		}
		pcolStats[j]=pregFuns[regFun]((pcolArr + j*nrb),nrb);
	}
	for(int i = 0; i<nrb;i++){
		prowStats[i]=pregFuns[regFun]((prowArr + i*ncb),ncb);
	}
	free(prowArr);
	free(pcolArr);
	
	for(int j = 0; j<ncb; j++){
		for(int i = 0; i<nrb; i++){
			res+= max(preSpecVal - min(pcolStats[j], prowStats[i]),0.0);
		}
	}
	
	if(*pmulReg==1) {
		free(prowStats);
		free(pcolStats);		
		return(res);
	} else {
		double res2=0;
		for(int j = 0; j<ncb; j++){
			res2+= max(preSpecVal - pcolStats[j],0.0);
		}
		for(int i = 0; i<nrb; i++){
			res2+= max(preSpecVal - prowStats[i],0.0);
		}
		free(prowStats);
		free(pcolStats);

		return(min(res,res2));

	}
}


/* a function for computing error of  the row-dominant block - binary blockmodeling*/
double valRdo(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
	mrd - minimal row deviation
*/

	int baseInd=relN*nr*nc;
	double prs=0;
	double mrd=0;
	

	for(int i = 0; i<nrb; i++){
		prs=0;
		for(int j = 0; j<ncb; j++){
			prs += max(preSpecVal-pM[baseInd+nc*pColInd[j]+pRowInd[i]],0.0);
		}
		if(mrd>prs){mrd = prs; }
	}
	
	if(*pmulReg==1) {
		return((mrd)*ncb);
	} else {
		return((mrd));
	}
}


/* a function for computing error of  the row-dominant - binary blockmodeling - diagonal*/
double valRdoDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/
	int baseInd=relN*nr*nc;
	double prs=0;
	double mrd=0;
	double resDiag=0.0;
	double tmp;

	
	for(int i = 0; i<nrb; i++){
		resDiag += max(preSpecVal-pM[baseInd+ nc*pColInd[i]+pRowInd[i]],0.0);
	}
	for(int i = 0; i<nrb; i++){
		for(int j = 0; j<ncb; j++){	
			tmp = max(preSpecVal-pM[baseInd+nc*pColInd[j]+pRowInd[i]],0.0);
			if((i==j) && (tmp > resDiag)){
				prs += resDiag;
			}else{
				prs +=tmp;
			}
			
		}
		
		if(mrd>prs){mrd = prs;}
	}
	
	if(*pmulReg==1) {
		return((mrd)*ncb);
	} else {
		return((mrd));
	}
}

/* a function for computing error of  the row-dominant block - binary blockmodeling - ignore diagonal*/
double valRdoIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
	mrd - minimal row deviation
*/

	int baseInd=relN*nr*nc;
	double prs=0;
	double mrd=0;
	

	for(int i = 0; i<nrb; i++){
		prs=0;
		for(int j = 0; j<ncb; j++){
			if(i!=j){ 
				prs += max(preSpecVal-pM[baseInd+nc*pColInd[j]+pRowInd[i]],0.0);
			}
		}
		if(mrd>prs){mrd = prs; }
	}
	
	if(*pmulReg==1) {
		return((mrd)*ncb);
	} else {
		return((mrd));
	}
}



/* a function for computing error of  the row-dominant block - binary blockmodeling*/
double valCdo(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
	mcd - minimal row deviation
*/
	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	double mcd=0;

	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcs += max(preSpecVal-pM[ind2d+pRowInd[i]],0.0);
		}
		if(mcd<pcs){mcd = pcs; }
	}
	
	if(*pmulReg==1) {
		return((mcd)*nrb);
	} else {
		return((mcd));
	}
}


/* a function for computing error of  the row-dominant - binary blockmodeling - diagonal*/
double valCdoDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	
/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	double mcd=0;
	double resDiag=0.0;
	double tmp;

	for(int i = 0; i<nrb; i++){
		resDiag += max(preSpecVal-pM[baseInd+ nc*pColInd[i]+pRowInd[i]],0.0);
	}


	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			tmp = max(preSpecVal-pM[ind2d+pRowInd[i]],0.0);
			if((i==j) && (tmp > resDiag)){
				pcs += resDiag;
			}else{
				pcs +=tmp;
			}
		}
		if(mcd<pcs){mcd = pcs; }
	}
	
	if(*pmulReg==1) {
		return((mcd)*nrb);
	} else {
		return((mcd));
	}
}

/* a function for computing error of  the row-dominant block - binary blockmodeling - ignore diagonal*/
double valCdoIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

/*
	nr - number of rows in the whole matrix/network
	nc - number of columns in the whole matrix/network
	nrb - number of rows in the block
	ncb - number of columns in the block
	mcd - minimal row deviation
*/

	int baseInd=relN*nr*nc;
	int ind2d;
	double pcs=0;
	double mcd=0;

	for(int j = 0; j<ncb; j++){
		pcs=0;
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			if(i!=j){ 
				pcs += max(preSpecVal-pM[ind2d+pRowInd[i]],0.0);
			}
		}
		if(mcd<pcs){mcd = pcs; }
	}
	
	if(*pmulReg==1) {
		return((mcd)*nrb);
	} else {
		return((mcd));
	}

	
}




/* a function for computing error of  the column-regular block - valued blockmodeling*/
double valCre(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *pcolArr;
	pcolArr = (double *) malloc(ncb*nrb*sizeof(double));

	double colStats;


	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcolArr[j*nrb + i] = pM[ind2d+pRowInd[i]];
		}
		colStats=pregFuns[regFun]((pcolArr + j*nrb),nrb);
		
		if(*pmulReg==1) {
			res+= max(preSpecVal - colStats,0.0)*nrb;
		} else {
			res+= max(preSpecVal - colStats,0.0);
		}
	}
	free(pcolArr);

	return(res);
}










/* a function for computing error of  the row-regular block - valued blockmodeling*/
double valRre(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *prowArr;
	prowArr = (double *) malloc(nrb*ncb*sizeof(double));

	double rowStats;



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			prowArr[i*ncb + j] = pM[ind2d+pRowInd[i]];
		}
	}
	for(int i = 0; i<nrb;i++){
		rowStats=pregFuns[regFun]((prowArr + i*ncb),ncb);
		
		if(*pmulReg==1) {
			res+= max(preSpecVal -  rowStats,0.0)*ncb;
		} else {
			res+= max(preSpecVal -  rowStats,0.0);
		}
		
	}
	free(prowArr);

	
	return(res);
}

/* a function for computing error of  the column-functional block - valued blockmodeling*/
double valCfn(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *pcolArr;
	pcolArr = (double *) malloc(ncb*nrb*sizeof(double));

	double colStats;
	double colSums;



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			pcolArr[j*nrb + i] = pM[ind2d+pRowInd[i]];
		}
		colStats=maxv((pcolArr + j*nrb),nrb);
		colSums=sumv((pcolArr + j*nrb),nrb);
		if(*pmulReg==1) {
			res+= max(preSpecVal - colStats,0.0)*nrb + colSums -colStats;
		} else {
			res+= max(preSpecVal - colStats,0.0) + colSums -colStats;
		}		
	}
	free(pcolArr);
	return(res);
}





/* a function for computing error of  the row-functional block - valued blockmodeling*/
double valRfn(const double *pM, const int nr, const int nc, const int relN, const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){
	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	double *prowArr;
	prowArr = (double *) malloc(nrb*ncb*sizeof(double));

	double rowStats;
	double rowSums;



	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			prowArr[i*ncb + j] = pM[ind2d+pRowInd[i]];
		}
	}
	for(int i = 0; i<nrb;i++){
		rowStats=maxv((prowArr + i*ncb),ncb);
		rowSums=sumv((prowArr + i*ncb),ncb);

		
		if(*pmulReg==1) {
			res+= max(preSpecVal -  rowStats,0.0)*ncb + rowSums-rowStats; 
		} else {
			res+= max(preSpecVal -  rowStats,0.0) + rowSums-rowStats; 
		}	
	}
	free(prowArr);
	
	return(res);
}







/* a function for computing error of  the average/density block - values/binary blockmodeling*/
double valAvg(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0.0;
	int baseInd=relN*nr*nc;
	int ind2d;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			res += pM[ind2d+pRowInd[i]];
		}
	}
	return(max(0,preSpecVal*ncb*nrb - res));
}

/* a function for computing error of  the average/density block - values/binary blockmodeling - diagonal*/
double valAvgDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	double diagRes=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	if(nrb==1){
		return(0.0);
	} else {
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			diagRes += pM[ind2d+pRowInd[j]];
			for(int i = (j + 1); i<nrb; i++){
				res += pM[ind2d+pRowInd[i]];
				res += pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
			}
		}
		return(max(0,preSpecVal*ncb*(nrb - 1) - res) + min(diagRes,(preSpecVal*nrb-diagRes)));
	}
}


/* a function for computing error of  the average/density block - values/binary blockmodeling -  diagonal ignore*/
double valAvgIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	double res=0;
	int baseInd=relN*nr*nc;
	int ind2d;

	if(nrb==1){
		return(0);
	} else {
		for(int j = 0; j<ncb; j++){
			ind2d=baseInd+ nc*pColInd[j];
			for(int i = (j + 1); i<nrb; i++){
				res += pM[ind2d+pRowInd[i]];
				res += pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
			}
		}
		return(max(0,preSpecVal*ncb*(nrb - 1) - res));
	}
}


/* a function for computing error of  the complete block - valued blockmodeling*/
double valCom(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	int baseInd=relN*nr*nc;
	int ind2d;
	double res=0.0;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			res+=max(preSpecVal-pM[ind2d+pRowInd[i]],0.0);
		}
	}
	return(res);
}

/* a function for computing error of  the complete block - valued blockmodeling - diagonal*/
double valComDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	int baseInd=relN*nr*nc;
	int ind2d;
	double res=0.0;
	double resDiag=0.0;
	double resDiag0=0.0;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];

		resDiag += max(preSpecVal-pM[ind2d+pRowInd[j]],0.0);
		resDiag0 += pM[ind2d+pRowInd[j]];

		for(int i = (j + 1); i<nrb; i++){
			res+=max(preSpecVal-pM[ind2d+pRowInd[i]],0.0);
			res+=max(preSpecVal-pM[baseInd+ nc*pColInd[i] + pRowInd[j]],0.0);
		}

	}
	return(res+min(resDiag,resDiag0));
}


/* a function for computing error of  the complete block - valued blockmodeling - ignore diagonal*/
double valComIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	int baseInd=relN*nr*nc;
	int ind2d;
	double res=0.0;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];

		for(int i = (j + 1); i<nrb; i++){
			res+=max(preSpecVal-pM[ind2d+pRowInd[i]],0.0);
			res+=max(preSpecVal-pM[baseInd+ nc*pColInd[i] + pRowInd[j]],0.0);
		}

	}
	return(res);
}

/* a function for computing error of  the null block - valued blockmodeling*/
double valNul(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	int baseInd=relN*nr*nc;
	int ind2d;
	double res=0.0;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = 0; i<nrb; i++){
			res+=pM[ind2d+pRowInd[i]];
		}
	}
	return(res);
}

/* a function for computing error of  the null block - valued blockmodeling - diagonal*/
double valNulDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	int baseInd=relN*nr*nc;
	int ind2d;
	double res=0.0;
	double resDiag=0.0;
	double resDiag0=0.0;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];

		resDiag += max(preSpecVal-pM[ind2d+pRowInd[j]],0.0);
		resDiag0 += pM[ind2d+pRowInd[j]];

		for(int i = (j + 1); i<nrb; i++){
			res+= pM[ind2d+pRowInd[i]];
			res+= pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
		}

	}
	return(res+min(resDiag,resDiag0));
}


/* a function for computing error of  the null block - valued blockmodeling - ignore diagonal*/
double valNulIgnoreDiag(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal, const double preSpecVal, const int *pmulReg){

	int baseInd=relN*nr*nc;
	int ind2d;
	double res=0.0;

	for(int j = 0; j<ncb; j++){
		ind2d=baseInd+ nc*pColInd[j];
		for(int i = (j + 1); i<nrb; i++){
			res+= pM[ind2d+pRowInd[i]];
			res+= pM[baseInd+ nc*pColInd[i] + pRowInd[j]];
		}

	}
	return(res);
}


/* Definition of an array of pointers to a function for computing block errors*/
double (*pBlockErr[nApproaches][nBlockTypes][3])(const double *pM, const int nr, const int nc, const int relN,const int nrb,const int ncb,const int *pRowInd, const int *pColInd, const int regFun, const int homFun, const int usePreSpecVal,const double preSpecVal, const int *pmulReg);




void critFun(const double *pM, const int *pnr, const int *pnc,  const int *pnRel, const int *pisTwoMode, const int *pisSym,const int *pdiag, const int *pnColClus, const int *pnRowClus, const int *pnUnitsRowClu, const int *pnUnitsColClu, const int *prowParArr, const int *pcolParArr,const int *papproaches, const int *pmaxBlockTypes,const int *pnBlockTypeByBlock, const int *pblocks, int *pIM, double *pEM, double *pEarr, double *perr, const int *pjustChange, const int *prowCluChange, const int *pcolCluChange, const int *psameIM, const int *pregFun, const int *phomFun, const int *pusePreSpec, const double *ppreSpecM, const double *pcombWeights, const int *pmulReg, const int *pnrInSetByClusters){
/*
double *pM - pointer to array or matrix representing the (multirotational) network
int *pnr - pointer to the number of rows
int *pnc - pointer to the number of columns
int *pisTwoMode - pointer to 0 (false) or 1 (true) specifying it the network is two-mode
int *pisSym - pointer to array of length (nRel - number of relation) specifying if the matrix (for each relation) is symetric) (0 - as any other value, 1 - seperately, 2 - ignore)
int *pdiag - pointer to array of length (nRel - number of relation) 0 (diag the same), 1 (diag special) or 2 (ignore values on diag) specifying how to treat the diagonal elements
int *pnRel - pointer to the number of relations
int *pnColClus - pointer to the number of column clusters
int *pnRowClus - pointer to the number of column clusters
int *pnUnitsRowClu - pointer to the array of the number of members of each row cluster
int *pnUnitsColClu - pointer to the array of the number of members of each row cluster
int *prowParArr - pointer to the array of arrays (one for each row cluster) of members of each row cluster
int *pcolParArr - pointer to the array of arrays (one for each col cluster) of members of each col cluster
int *papproaches - pointer to the array specifying approach - one for each relation
int *pmaxBlockTypes - pointer to maximum number of used block types
int *pnBlockTypeByBlock - pointer to 3d array (Rel, row, col) specifying the number of used allowed block types
int *pblocks - pointer to the 4d array (pmaxBlockTypes, Rel, row, col) specifying allowed block types
### Not implemetned ### double *pblocks - pointer to error multipliers by blocks
int *pIM - pointer to 3d array (Rel, row, col) specifying the image matrix
double *pEM - pointer to 3d array (Rel, row, col) specifying the error for each block
double *pEarr - pointer to the 4d array (pmaxBlockTypes, Rel, row, col) specifying the errors for each allowed block type - it is very important that the value is Infinitive for block types that are not allowed
double *perr - pointer to the error returned by this function
int *pjustChange - pointer to a value specifying if only the errors for changed clusters should be computed
int *prowCluChange - pointer to an array holding the two row clusters where the change occurred
int *pcolCluChange - pointer to an array holding the col row clusters where the change occurred
int *psameIM - pointer to 0 (false) or 1 (true) specifying if the image has to be the same for all relations
int *pregFun - pointer to the 4d array (pmaxBlockTypes, Rel, row, col) specifying the "summary" function used in f-regular line blocks
int *phomFun - pointer to the array (one value for each rel) function used used for computing measure of variability in homogeneity blockmodeling
int *pusePreSpec - pointer to 4d array (pmaxBlockTypes, Rel, row, col) specifying weather a the pre-specified value should be used when computing inconsistency
double *ppreSpecM - pointer to 4d array (pmaxBlockTypes, Rel, row, col) specifying the pre-specified value to be used when computing inconsistency
double *pcombWeights - pointer to a array of weights of the same dimensions as blocks
const int *pmulReg - pointer to 0 (false) or 1 (true) specifying if the errors that apply to rows/columns (and not to cells) should be multiplied by number of rows/columns.
int *pnrInSetByClusters - pointer to the vector of sizes of sets for row clusters. The length of the vector is equal to the number of row clusters and gives the number of units in a set to which a given row cluster belongs to. Should be set to vector of 0s if ll for group sizes should not be computed. 
*/

/* a lot of arguments not used: *pisSym, *pisTwoMode,.... */
/* *pisTwoMode argument is probably/perhaps not needed. Its information could be previously incorporated into *pdiag argument */
	/*int N = (*pnr)*(*pnc);*/

/*Rprintf("critFun - start \n");*/

	int ind2d, ind3d, ind4d;
	double minBlockErrVal;
	int minBlockErrInd;
	int iRel=0, iDiag;
	int iBlockType=0;


	// Rprintf("pnrInSetByClusters: ");
	// for( int i=0;i<*pnRowClus;i++){
		// Rprintf("%i ", pnrInSetByClusters[i]);
	// }
	
	/* assiginig functions to the array of functions*/
	/* array of pointers to functions for f-regular (and similar) blocks
	usage: pregFuns[function]
	function: 	0 - maxv (max)
				1 - sumv (sum)
				2 - meanv (mean)
	 */

	pregFuns[0]=maxv;
	pregFuns[1]=sumv;
	pregFuns[2]=meanv;

/*Rprintf("critFun - 1\n");*/


	/* assiginig functions to the array of functions*/
	/* array of pointers to functions for computing measure of variability
	usage: phom[measureOfVariability][prespecifiedValue]
	measureOfVariability: 	0 - ss (sum of squared deviations from the prespecified value (default mean))
							1 - ad (absolute deviations from the prespecified value (default median))
							2 - bll (- binary log-likelihood)
	prespecifiedValue:	0 - default (mean or median, depending on the measureOfVariability)
						1 - as specified in ppreSpecM
						2 - 0
						3 - maximum of the value specified in ppreSpecM and the default (mean or median, depending on the measureOfVariability)
	 */

	phom[0][0] = ss;
	phom[0][1] = ssP;
	phom[0][2] = ss0;
	phom[0][3] = ssPmin;
	phom[1][0] = ad;
	phom[1][1] = adP;
	phom[1][2] = ad0;
	phom[1][3] = adPmin;
	phom[2][0] = bll;
	phom[2][1] = bllP;
	phom[2][2] = bll0;
	phom[2][3] = bllPmin;

/*Rprintf("critFun - 2\n");*/

	/* assigning functions to the array of functions*/
	/* array of pointers to functions for computing block errors
	usage: pBlockErr[blockmodelingApproach][blockType][treatDiagonal]
	blockmodelingApproach: 	0 - homgeneity blockmodeling
							1 - binary blockmodeling
							2 - valued blockmodeling
	blockType:	0 - null
				1 - complete
				2 - column-functional
				3 - row-functional
				4 - regular
				5 - column-regular
				6 - row-regular
				7 - average
				8 - do not care
	treatDiagonal:  0 - as any other value
	                1 - separately
	                2 - ignore diagonal
	 */

	pBlockErr[0][0][0]=homNul;
	pBlockErr[0][0][1]=homNulDiag;
	pBlockErr[0][0][2]=homNulIgnoreDiag;

	pBlockErr[0][1][0]=homCom;
	pBlockErr[0][1][1]=homComDiag;
	pBlockErr[0][1][2]=homComIgnoreDiag;

	pBlockErr[0][2][0]=infBlock;
	pBlockErr[0][2][1]=infBlock;
	pBlockErr[0][2][2]=infBlock;

	pBlockErr[0][3][0]=infBlock;
	pBlockErr[0][3][1]=infBlock;
	pBlockErr[0][3][2]=infBlock;
	
	pBlockErr[0][4][0]=homCfn;
	pBlockErr[0][4][1]=homCfn;
	pBlockErr[0][4][2]=homCfn;

	pBlockErr[0][5][0]=homRfn;
	pBlockErr[0][5][1]=homRfn;
	pBlockErr[0][5][2]=homRfn;
	
	pBlockErr[0][6][0]=homReg;
	pBlockErr[0][6][1]=homReg;
	pBlockErr[0][6][2]=homReg;

	pBlockErr[0][7][0]=homCre;
	pBlockErr[0][7][1]=homCre;
	pBlockErr[0][7][2]=homCre;

	pBlockErr[0][8][0]=homRre;
	pBlockErr[0][8][1]=homRre;
	pBlockErr[0][8][2]=homRre;

/*There is no difference between complete and "average" block for homogeneity blockmodeling*/	
	pBlockErr[0][9][0]=homCom;
	pBlockErr[0][9][1]=homComDiag;
	pBlockErr[0][9][2]=homComIgnoreDiag;	

	pBlockErr[0][10][0]=doNotCare;
	pBlockErr[0][10][1]=doNotCare;
	pBlockErr[0][10][2]=doNotCare;

	pBlockErr[1][0][0]=binNul;
	pBlockErr[1][0][1]=binNulDiag;
	pBlockErr[1][0][2]=binNulIgnoreDiag;

	pBlockErr[1][1][0]=binCom;
	pBlockErr[1][1][1]=binComDiag;
	pBlockErr[1][1][2]=binComIgnoreDiag;

	pBlockErr[1][2][0]=binCdo;
	pBlockErr[1][2][1]=binCdoDiag;
	pBlockErr[1][2][2]=binCdoIgnoreDiag;

	pBlockErr[1][3][0]=binRdo;
	pBlockErr[1][3][1]=binRdoDiag;
	pBlockErr[1][3][2]=binRdoIgnoreDiag;


	pBlockErr[1][4][0]=binCfn;
	pBlockErr[1][4][1]=binCfn;
	pBlockErr[1][4][2]=binCfn;

	pBlockErr[1][5][0]=binRfn;
	pBlockErr[1][5][1]=binRfn;
	pBlockErr[1][5][2]=binRfn;

	pBlockErr[1][6][0]=binReg;
	pBlockErr[1][6][1]=binReg;
	pBlockErr[1][6][2]=binReg;

	pBlockErr[1][7][0]=binCre;
	pBlockErr[1][7][1]=binCre;
	pBlockErr[1][7][2]=binCre;

	pBlockErr[1][8][0]=binRre;
	pBlockErr[1][8][1]=binRre;
	pBlockErr[1][8][2]=binRre;

	
/*Functions for density (binary) block and for average valued blocks are the same*/	
	pBlockErr[1][9][0]=valAvg;
	pBlockErr[1][9][1]=valAvgDiag;
	pBlockErr[1][9][2]=valAvgIgnoreDiag;
	
	pBlockErr[1][10][0]=doNotCare;
	pBlockErr[1][10][1]=doNotCare;
	pBlockErr[1][10][2]=doNotCare;

	pBlockErr[2][0][0]=valNul;
	pBlockErr[2][0][1]=valNulDiag;
	pBlockErr[2][0][2]=valNulIgnoreDiag;

	pBlockErr[2][1][0]=valCom;
	pBlockErr[2][1][1]=valComDiag;
	pBlockErr[2][1][2]=valComIgnoreDiag;
	
	pBlockErr[2][2][0]=valCdo;
	pBlockErr[2][2][1]=valCdoDiag;
	pBlockErr[2][2][2]=valCdoIgnoreDiag;

	pBlockErr[2][3][0]=valRdo;
	pBlockErr[2][3][1]=valRdoDiag;
	pBlockErr[2][3][2]=valRdoIgnoreDiag;	

	pBlockErr[2][4][0]=valCfn;
	pBlockErr[2][4][1]=valCfn;
	pBlockErr[2][4][2]=valCfn;

	pBlockErr[2][5][0]=valRfn;
	pBlockErr[2][5][1]=valRfn;
	pBlockErr[2][5][2]=valRfn;

	pBlockErr[2][6][0]=valReg;
	pBlockErr[2][6][1]=valReg;
	pBlockErr[2][6][2]=valReg;

	pBlockErr[2][7][0]=valCre;
	pBlockErr[2][7][1]=valCre;
	pBlockErr[2][7][2]=valCre;

	pBlockErr[2][8][0]=valRre;
	pBlockErr[2][8][1]=valRre;
	pBlockErr[2][8][2]=valRre;

	pBlockErr[2][9][0]=valAvg;
	pBlockErr[2][9][1]=valAvgDiag;
	pBlockErr[2][9][2]=valAvgIgnoreDiag;	

	pBlockErr[2][10][0]=doNotCare;
	pBlockErr[2][10][1]=doNotCare;
	pBlockErr[2][10][2]=doNotCare;



/*Rprintf("critFun - 3\n");*/

	double *pEMarrAllRel;
	pEMarrAllRel = 	(double *) malloc((*pmaxBlockTypes)*sizeof(double));

	if(*pjustChange){
/*Rprintf("JustChange=%i\n", *pjustChange);*/
		/* is it justified to have two options for that */
		*perr=0.0;
		for(int iColClu=0;iColClu<*pnColClus;iColClu++){	
			int colChange = (iColClu==pcolCluChange[0])|(iColClu==pcolCluChange[1]);
			for(int iRowClu=0;iRowClu<*pnRowClus;iRowClu++){
				if(colChange | (iRowClu==prowCluChange[0])|(iRowClu==prowCluChange[1])){
					ind2d=iColClu*(*pnRowClus) + iRowClu;

					for(iRel=0; iRel<(*pnRel);iRel++){
						iDiag = (iColClu==iRowClu)? pdiag[iRel]: 0;
						ind3d= (ind2d*(*pnRel)+ iRel);
						minBlockErrVal = INFINITY;
						minBlockErrInd=0;
						for(iBlockType=0;iBlockType<(pnBlockTypeByBlock[ind3d]);iBlockType++){

							ind4d=ind3d*(*pmaxBlockTypes)+iBlockType;
							pEarr[ind4d]=pcombWeights[ind4d]*pBlockErr[papproaches[iRel]][pblocks[ind4d]][iDiag](pM,*pnr,*pnc,iRel,pnUnitsRowClu[iRowClu],pnUnitsColClu[iColClu],(prowParArr + iRowClu*(*pnr)),(pcolParArr +iColClu*(*pnc)),pregFun[ind4d],phomFun[iRel],pusePreSpec[ind4d],ppreSpecM[ind4d], pmulReg);

							if((!(*psameIM))&& (pEarr[ind4d]<minBlockErrVal)){
								minBlockErrVal=pEarr[ind4d];
								minBlockErrInd=ind4d;
							}
						}
						if(!(*psameIM)){
							pEM[ind3d]=minBlockErrVal;
							pIM[ind3d]=pblocks[minBlockErrInd];
							*perr+=minBlockErrVal;
						}
					}
					if(*psameIM){
						minBlockErrVal=INFINITY;
						minBlockErrInd=0;
						for(iBlockType=0;iBlockType<(*pmaxBlockTypes);iBlockType++){
							pEMarrAllRel[iBlockType]=0.0;
							for(int iRel=0; iRel<(*pnRel);iRel++){
								pEMarrAllRel[iBlockType]+=pEarr[(ind2d*(*pnColClus)+ iRel)*(*pmaxBlockTypes)+iBlockType];
							}
							if(pEMarrAllRel[iBlockType]<minBlockErrVal){
								minBlockErrVal=pEMarrAllRel[iBlockType];
								minBlockErrInd=iBlockType;
							}
						}

						for(iRel=0; iRel<(*pnRel);iRel++){
							ind3d= (ind2d*(*pnColClus)+ iRel);
							pEM[ind3d]=pEarr[ind3d + minBlockErrInd];
							pIM[ind3d]=pblocks[ind3d + minBlockErrInd];
							*perr+=pEM[ind3d];
						}
					}
				} else {
					ind2d=iColClu*(*pnRowClus) + iRowClu;
					for(iRel=0; iRel<(*pnRel);iRel++){
						ind3d= (ind2d*(*pnRel)+ iRel);
						*perr+=pEM[ind3d];
					}
				}
			}
		}
	}else{
/*Rprintf("critFun - 3.5\n");*/
		*perr=0.0;
		for(int iColClu=0;iColClu<*pnColClus;iColClu++){
			for(int iRowClu=0;iRowClu<*pnRowClus;iRowClu++){
				ind2d=iColClu*(*pnRowClus) + iRowClu;

				for(iRel=0; iRel<(*pnRel);iRel++){
					iDiag = (iColClu==iRowClu)? pdiag[iRel]: 0;
					ind3d= (ind2d*(*pnRel)+ iRel);
					minBlockErrVal = INFINITY;
					minBlockErrInd=0;
					for(iBlockType=0;iBlockType<(pnBlockTypeByBlock[ind3d]);iBlockType++){
						ind4d=ind3d*(*pmaxBlockTypes)+iBlockType;
/*Rprintf("critFun - 4\n");*/
/*Rprintf("iColClu = %i, iRowClu = %i, iRel=%i\n", iColClu, iRowClu, iRel);*/
/*Rprintf("approach = %i, blockType = %i, iDiag = %i\n", papproaches[iRel], pblocks[ind4d], iDiag);*/
/*Rprintf("ind4d = %i\n", ind4d);*/
						/* double temp */
						pEarr[ind4d]=pcombWeights[ind4d] * pBlockErr[papproaches[iRel]][pblocks[ind4d]][iDiag](pM,*pnr,*pnc,iRel,pnUnitsRowClu[iRowClu],pnUnitsColClu[iColClu],(prowParArr + iRowClu*(*pnr)),(pcolParArr +iColClu*(*pnc)),pregFun[ind4d],phomFun[iRel],pusePreSpec[ind4d],ppreSpecM[ind4d], pmulReg);
/*Rprintf("blockErr = %.2f\n", temp);*/
/*Rprintf("critFun - 4.9\n");*/
						/* pEarr[ind4d]=temp; */
/*Rprintf("critFun - 5\n");*/
						if((!(*psameIM))&& (pEarr[ind4d]<minBlockErrVal)){
							minBlockErrVal=pEarr[ind4d];
							minBlockErrInd=ind4d;
						}
					}
					if(!(*psameIM)){
						pEM[ind3d]=minBlockErrVal;
						pIM[ind3d]=pblocks[minBlockErrInd];
						*perr+=minBlockErrVal;
					}
				}
				if(*psameIM){
					minBlockErrVal=INFINITY;
					minBlockErrInd=0;
					for(iBlockType=0;iBlockType<(*pmaxBlockTypes);iBlockType++){
						pEMarrAllRel[iBlockType]=0.0;
						for(int iRel=0; iRel<(*pnRel);iRel++){
							pEMarrAllRel[iBlockType]+=pEarr[(ind2d*(*pnColClus)+ iRel)*(*pmaxBlockTypes)+iBlockType];
						}
						if(pEMarrAllRel[iBlockType]<minBlockErrVal){
							minBlockErrVal=pEMarrAllRel[iBlockType];
							minBlockErrInd=iBlockType;
						}
					}

					for(iRel=0; iRel<(*pnRel);iRel++){
						ind3d= (ind2d*(*pnColClus)+ iRel);
						pEM[ind3d]=pEarr[ind3d + minBlockErrInd];
						pIM[ind3d]=pblocks[ind3d + minBlockErrInd];
						*perr+=pEM[ind3d];
					}
				}
			}
		}
	}
	
// Rprintf("err = %.3f\n", *perr);
	// Rprintf("pnrInSetByClusters: ");
	// for( int i=0;i<*pnRowClus;i++){
		// Rprintf("%i ", pnrInSetByClusters[i]);
	// }
	// Rprintf("\npnUnitsRowClu: ");
	// for( int i=0;i<*pnRowClus;i++){
		// Rprintf("%i ", pnUnitsRowClu[i]);
	// }	
	// Rprintf("\npnrInSetByClusters[0] = %i\n", pnrInSetByClusters[0]);
	if(pnrInSetByClusters[0]>0){
	// Rprintf("True\n");
		// for(int iColClu=0;iColClu<*pnColClus;iColClu++){
		for(int iRowClu=0;iRowClu<*pnRowClus;iRowClu++){
// Rprintf("iRowClu = %i\n", iRowClu);
// Rprintf("pnUnitsRowClu[iRowClu] = %i\n", pnUnitsRowClu[iRowClu]);
// Rprintf("pnrInSetByClusters[iRowClu] = %i\n", pnrInSetByClusters[iRowClu]);		
			*perr += -pnUnitsRowClu[iRowClu]*log(pnUnitsRowClu[iRowClu]*1.0/pnrInSetByClusters[iRowClu]);
			// Rprintf("err = %.3f\n", *perr);
		}
	}
	free(pEMarrAllRel);
/*Rprintf("critFun - end \n");*/
}



/* the function below converts an array representation of a partition to a vector representation of a partition */
void parArr2Vec(const int *pn, const int *pnClus, const int *pnUnitsClu, const int *pParArr, int *pParVec){
	/*pParVec = (int *) malloc((*pn)*sizeof(int));*/
	for(int iClu=0;iClu<*pnClus;iClu++){
		for(int iCluUnit=0;iCluUnit<pnUnitsClu[iClu];iCluUnit++){
			pParVec[pParArr[iClu*(*pn)+iCluUnit]]=iClu;
		}
	}
}

/* the function below converts a vector representation of a partition to an array representation of a partition */
void parVec2Arr(const int *pn, int *pnClus, int *pnUnitsClu, int *pParArr, const int *pParVec){
/*	Rprintf("OK1"); */
	int nClus=0;
	for(int i=0;i<*pn;i++){
		if(pParVec[i]>=nClus) nClus = pParVec[i]+1;
	}
/*	Rprintf("OK2"); */
	*pnClus = nClus;
/*	Rprintf("OK3"); */
	/*pnUnitsClu = (int *) malloc((*pnClus)*sizeof(int));*/
	/*pParArr = (int *) malloc((*pnClus)*(*pn)*sizeof(int));*/
	for(int i=0;i<*pn;i++){
		pParArr[pParVec[i]*(*pn)+pnUnitsClu[pParVec[i]]]=i;
		pnUnitsClu[pParVec[i]]++;
		Rprintf("OK4.%i", i);
	}
/*	Rprintf("OK5"); */
}




void updateResults(const int *pnc, const int *pnRel, const int *pnColClus, const int *pnRowClus, const int *pmaxBlockTypes,  const int *psourcenUnitsRowClu, const int *psourcerowParArr, const int *psourceIM, const double *psourceEM, const double *psourceEarr, const double *psourceerr, int *pdestnUnitsRowClu, int *pdestrowParArr, int *pdestIM, double *pdestEM, double *pdestEarr, double *pdesterr){
	/*update dest results */

	*pdesterr = *psourceerr;

	for(int i=0;i<*pnRowClus;i++){
		pdestnUnitsRowClu[i] = psourcenUnitsRowClu[i];
	}
	for(int i=0;i<((*pnRowClus)*(*pnc));i++){
		pdestrowParArr[i] = psourcerowParArr[i];
	}

	/* image matrix */
	for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
		pdestIM[i] = psourceIM[i];
	}

	/* number of block types by block - not needed
	int *pdestnBlockTypeByBlock;
	pdestnBlockTypeByBlock = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
	for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
		pdestnBlockTypeByBlock[i] = pnBlockTypeByBlock[i];
	} */


	/* error matrix */
	for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
		pdestEM[i] = psourceEM[i];
	}

	/* error array by block types*/
	for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
		pdestEarr[i] = psourceEarr[i];
	}
}

/* for now this function moves to improved partition as soon as it finds one */
/* however, the "move" is selected randomly, while it is true that "moves" are tried before "exchanges" */
void optPar(const double *pM, const int *pnr, const int *pnc,  const int *pnRel, const int *pisTwoMode, const int *pisSym, const int *pdiag, const int *pjustMove, const int *pnImp, const int *pnColClus, const int *pnRowClus, int *pnUnitsRowClu, int *pnUnitsColClu, int *prowParArr, int *pcolParArr,const int *papproaches, const int *pmaxBlockTypes,const int *pnBlockTypeByBlock, const int *pblocks, int *pIM, double *pEM, double *pEarr, double *perr, const int *pjustChange, int *prowCluChange, int *pcolCluChange, const int *psameIM, const int *pregFun, const int *phomFun, const int *pusePreSpec, const double *ppreSpecM, const int *pminUnitsRowCluster, const int *pminUnitsColCluster, const int *pmaxUnitsRowCluster, const int *pmaxUnitsColCluster, int *psameErr, int *pnIter, const double *pcombWeights, const int *pexchageClusters, const int *pmulReg, const int *pnrInSetByClusters){
	/*
	double *pM - pointer to array or matrix representiing the (multirelational) network
	int *pnr - pointer to the number of rows
	int *pnc - pointer to the number of columns
	int *pisTwoMode - pointer to 0 (false) or 1 (true) specifying it the network is two-mode
	int *pisSym - pointer to array of length (nRel - number of relation) specifying if the matrix (for each relation) is symmetric) 0 (false) or 1 (true) 
	int *pdiag - pointer to array of length (nRel - number of relation) specifying how to treat the diagonal elements (0 - as any other value, 1 - separately, 2 - ignore)
	int *pjustMove - pointer weather to try just moves (1) or moves and exchanges (0)
	int *pnImp - pointer to the number of improvements that must be found before before we move to a better partition.
	int *pnRel - pointer to the number of relations
	int *pnColClus - pointer to the number of column clusters
	int *pnRowClus - pointer to the number of column clusters
	int *pnUnitsRowClu - pointer to the array of the number of members of each row cluster
	int *pnUnitsColClu - pointer to the array of the number of members of each col cluster
	int *prowParArr - pointer to the array of arrays (one for each row cluster) of members of each row cluster
	int *pcolParArr - pointer to the array of arrays (one for each col cluster) of members of each col cluster
	int *papproaches - pointer to the array specifying approach - one for each relation
	int *pmaxBlockTypes - pointer to maximum number of used block types
	int *pnBlockTypeByBlock - pointer to 3d array (Rel, row, col) specifying the number of used allowed block types
	int *pblocks - pointer to the 4d array (nBlockTypesByBlock, Rel, row, col) specifying allowed block types
	int *pIM - pointer to 3d array (Rel, row, col) specifying the image matrix
	double *pEM - pointer to 3d array (Rel, row, col) specifying the error for each block
	double *pEarr - pointer to the 4d array ((*pmaxBlockTypes), Rel, row, col) specifying the errors for each allowed block type - it is very important that the value is Infinitive for block types that are not allowed
	double *perr - pointer to the total error
	int *pjustChange - pointer to a value specifying if only the errors for changed clusters should be computed
	int *prowCluChange - pointer to an array holding the two row clusters where the change occurred
	int *pcolCluChange - pointer to an array holding the col row clusters where the change occurred
	int *psameIM - pointer to 0 (false) or 1 (true) specifying if the image has to be the same for all relations
	int *pregFun - pointer to the 4d array ((*pmaxBlockTypes), Rel, row, col) specifying the "summary" function used in f-regular line blocks
	int *phomFun - pointer to the array (one value for each rel) function used used for computing measure of variability in sum of squares blockmodeling
	int *pusePreSpec - pointer to 4d array ((*pmaxBlockTypes), Rel, row, col) specifying weather a the pre-specified value should be used when computing inconsistency
	double *ppreSpecM - pointer to 4d array ((*pmaxBlockTypes), Rel, row, col) specifying the pre-specified value to be used when computing inconsistency
	int *pminUnitsRowCluster - pointer to the minimum number of units in row cluster
	int *pminUnitsColCluster - pointer to the minimum number of units in col cluster
	int *pmaxUnitsRowCluster - pointer to the maximum number of units in row cluster
	int *pmaxUnitsColCluster - pointer to the maximum number of units in col cluster
	double *pcombWeights - pointer to a array of weights of the same dimensions as blocks
	int *pexchageClusters - pointer to a matrix (nRowClust, nColClus) showing which clusters are exchangeable
	int *pnrInSetByClusters - pointer to the vector of sizes of sets for row clusters. The length of the vector is equal to the number of row clusters and gives the number of units in a set to which a given row cluster belongs to. Should be set to vector of 0s if ll for group sizes should not be computed. 
	*/

/*Rprintf("OptParC\n");*/
/**/
/*Rprintf("pM: ");*/
/*for( int i=0;i<(*pnr)*(*pnc)*(*pnRel);i++){*/
/*	Rprintf("%f ", pM[i]);*/
/*}*/
/*Rprintf("\n");*/

/*	int *pzero;
	pzero = (int *) malloc(sizeof(int));
	*pzero = 0; */
	int zero = 0;

/*	Rprintf("test1");*/
    GetRNGstate();	/* Get .Random.seed from R */
	if(*pisTwoMode){
		Rprintf("Optimization of two-mode networks is not yet supported\n");
	} else {

		critFun(pM, pnr, pnc,  pnRel, pisTwoMode, pisSym,  pdiag, pnColClus, pnRowClus, pnUnitsRowClu, pnUnitsColClu, prowParArr, pcolParArr, papproaches, pmaxBlockTypes, pnBlockTypeByBlock, pblocks, pIM, pEM, pEarr, perr, &zero, prowCluChange, pcolCluChange, psameIM, pregFun, phomFun, pusePreSpec,  ppreSpecM, pcombWeights, pmulReg, pnrInSetByClusters);
/*Rprintf("Initial error = %.2f\n", *perr);*/

		/* prepare temporary objects - start*/


		/* best result  - start*/

		/* partition*/
		int *pbestrowParArr;
		int *pbestnUnitsRowClu;
		pbestnUnitsRowClu = (int *) malloc((*pnRowClus)*sizeof(int));
		pbestrowParArr = (int *) malloc((*pnRowClus)*(*pnc)*sizeof(int));
		for(int i=0;i<*pnRowClus;i++){
			pbestnUnitsRowClu[i] = pnUnitsRowClu[i];
		}
		for(int i=0;i<((*pnRowClus)*(*pnc));i++){
			pbestrowParArr[i] = prowParArr[i];
		}

		/* image matrix */
		int *pbestIM;
		pbestIM = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestIM[i] = pIM[i];
		}

		/* number of block types by block - not needed
		int *pbestnBlockTypeByBlock;
		pbestnBlockTypeByBlock = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestnBlockTypeByBlock[i] = pnBlockTypeByBlock[i];
		} */


		/* error matrix */
		double *pbestEM;
		pbestEM = (double *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestEM[i] = pEM[i];
		}

		/* error array by block types*/
		double *pbestEarr;
		pbestEarr = (double *) malloc((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestEarr[i] = pEarr[i];
		}


		double *pbesterr;
		pbesterr = (double *) malloc(sizeof(double));
		*pbesterr = *perr;

		/* best result  - end*/




		/* temp result  - start*/

		/* partition*/
		int *ptemprowParArr;
		int *ptempnUnitsRowClu;
		ptempnUnitsRowClu = (int *) malloc((*pnRowClus)*sizeof(int));
		ptemprowParArr = (int *) malloc((*pnRowClus)*(*pnc)*sizeof(int));
		for(int i=0;i<*pnRowClus;i++){
			ptempnUnitsRowClu[i] = pnUnitsRowClu[i];
		}
		for(int i=0;i<((*pnRowClus)*(*pnc));i++){
			ptemprowParArr[i] = prowParArr[i];
		}


		/* image matrix */
		int *ptempIM;
		ptempIM = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempIM[i] = pIM[i];
		}

		/* number of block types by block - not needed
		int *ptempnBlockTypeByBlock;
		ptempnBlockTypeByBlock = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempnBlockTypeByBlock[i] = pnBlockTypeByBlock[i];
		} */


		/* error matrix */
		double *ptempEM;
		ptempEM = (double *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempEM[i] = pEM[i];
		}

		/* error array by block types*/
		double *ptempEarr;
		ptempEarr = (double *) malloc((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempEarr[i] = pEarr[i];
		}

		double *ptemperr;
		ptemperr = (double *) malloc(sizeof(double));
		*ptemperr = *perr;

		/* temp result  - end*/



		/* prepare temporary objects - end*/


		*psameErr = 1;
		int improve=1;
/*Rprintf("OK1\n");*/

		/* loop until no improvement is found */
		*pnIter=0;
		while(improve){
			*pnIter = *pnIter + 1;
			/* copy temp results to permanent  - start*/
			/* partition*/
			for(int i=0;i<*pnRowClus;i++){
				pnUnitsRowClu[i] = pbestnUnitsRowClu[i];
			}
			for(int i=0;i<((*pnRowClus)*(*pnc));i++){
				prowParArr[i] = pbestrowParArr[i];
			}

			/* image matrix */
			for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
				pIM[i] = pbestIM[i];
			}

			/* error matrix */
			for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
				pEM[i] = pbestEM[i];
			}

			/* error array by block types*/
			for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
				pEarr[i] = pbestEarr[i];
			}

			*perr = *pbesterr;

			/* copy temp results to permanent  - end*/

			improve=0;

			/*update temp results */
			updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr);


			/* to make the order of evaluation random - start */
			/*  randomize(); does not work */
			int rnd;
			int rndClusters[*pnRowClus];
			for(int i=0;i<*pnRowClus;i++){
				rndClusters[i]=i;
			}
			/* to make the order of evaluation random - end */

			int iClu, iClu2, iUnit, iUnit2;

			/* a loop over all clusters - random order */
			for(int iRndClu=0;iRndClu<*pnRowClus;iRndClu++){
/* Rprintf("Start loop cluster 1\n"); */

				/* to make the order of evaluation random - start*/
				rnd=randomInt(*pnRowClus-iRndClu);
				iClu=rndClusters[rnd];
				prowCluChange[0]=iClu;

				rndClusters[rnd]=rndClusters[*pnRowClus-iRndClu-1];
				/* to make the order of evaluation random - end*/

				/* a loop over all units inside clusters*/

				/* to make the order of evaluation random - start*/
				int rndUnitsInClu[pnUnitsRowClu[iClu]];
				for(int i=0;i<pnUnitsRowClu[iClu];i++){
					rndUnitsInClu[i]=i;
				}
				/* to make the order of evaluation random - end*/

/*Rprintf("rndUnitsInClu: ");*/
/*for( int i=0;i<pnUnitsRowClu[iClu];i++){*/
/*	Rprintf("%i ", rndUnitsInClu[i]);*/
/*}*/
/*Rprintf("\n");*/

				for(int iRndUnit=0;iRndUnit < pnUnitsRowClu[iClu];iRndUnit++){
/*Rprintf("Start loop unit in cluster 1\n");*/
					/* to make the order of evaluation random - start */
					rnd=randomInt(pnUnitsRowClu[iClu]-iRndUnit);
/*Rprintf("OK 1.01\n");*/
					iUnit=rndUnitsInClu[rnd];
/*Rprintf("OK 1.02\n");*/
					rndUnitsInClu[rnd]=rndUnitsInClu[pnUnitsRowClu[iClu] - iRndUnit - 1];
/*Rprintf("OK 1.03\n");*/
/*Rprintf("rnd = %i, pnUnitsRowClu[iClu]-iRndUnit= %i, iUnit=%i, rndUnitsInClu[rnd]=%i\n", rnd, pnUnitsRowClu[iClu]-iRndUnit, iUnit, rndUnitsInClu[rnd]);			*/
/*Rprintf("Selected unit ID = %i\n", ptemprowParArr[iClu*(*pnr)+iUnit]);*/
/*Rprintf("rndUnitsInClu: ");*/
/*for( int i=0;i<pnUnitsRowClu[iClu];i++){*/
/*	Rprintf("%i ", rndUnitsInClu[i]);*/
/*}*/
/*Rprintf("\n");*/

					/* to make the order of evaluation random - end*/


					/* to make the order of evaluation random - start*/
					int rndClusters2[*pnRowClus-1];
					for(int i=0;i<*pnRowClus;i++){
						if(i < iClu) {
							rndClusters2[i]=i;
						} else if (i>iClu){
							rndClusters2[i-1]=i;
						}
					}
/*Rprintf("OK 1.04\n");*/

					/* to make the order of evaluation random - end*/

					/* a loop over all other clusters - random order */
					for(int iRndClu2=0;iRndClu2<(*pnRowClus-1);iRndClu2++){
/*Rprintf("Start loop cluster 2\n");*/
						/* to make the order of evaluation random - start*/
						rnd=randomInt(*pnRowClus - 1 - iRndClu2);
/*Rprintf("rnd = %i, *pnRowClus - 1 - iRndClu2= %i\n", rnd, *pnRowClus - 1 - iRndClu2);			*/
						iClu2=rndClusters2[rnd];
						prowCluChange[1]=iClu2;
						rndClusters2[rnd]=rndClusters2[*pnRowClus - 2 - iRndClu2];
/*Rprintf("rndClusters2[rnd] = %i\n", rndClusters2[rnd]);*/
						if (!pexchageClusters[iClu*(*pnRowClus)+iClu2]){
							continue;
						}

						/* to make the order of evaluation random - end*/

						if((pnUnitsRowClu[iClu]>(*pminUnitsRowCluster))&&(pnUnitsRowClu[iClu2]<(*pmaxUnitsRowCluster))){
/*Rprintf("OK1.1\n");*/
							/* move unit to another cluster */
							ptemprowParArr[iClu2*(*pnr)+ptempnUnitsRowClu[iClu2]]=ptemprowParArr[iClu*(*pnr)+iUnit];
/*Rprintf("OK1.2\n");*/

							ptempnUnitsRowClu[iClu2]++;	/* this line must be after the above line */
/*Rprintf("OK1.3\n");*/
							ptempnUnitsRowClu[iClu]--;	/* this line must be before the line below */
/*Rprintf("OK1.4\n");*/
							ptemprowParArr[iClu*(*pnr)+iUnit]=ptemprowParArr[iClu*(*pnr)+ptempnUnitsRowClu[iClu]];

/*Rprintf("iClu = %i, iClu2= %i, iUnit=%i\n", iClu, iClu2, iUnit);*/
/*Rprintf("nClu = %i, nCluOld= %i, nClu2 = %i, nCluOld2= %i\n", ptempnUnitsRowClu[iClu], pnUnitsRowClu[iClu], ptempnUnitsRowClu[iClu2], pnUnitsRowClu[iClu2]); */
/*for(int i1=0;i1<(*pnRowClus);i1++){
	Rprintf("cluster = %i, unitsCluster= %i: ", i1, ptempnUnitsRowClu[i1]);
	for(int i2=0;i2<(ptempnUnitsRowClu[i1]);i2++){
		Rprintf("%i ", ptemprowParArr[i1*(*pnr)+i2]);
	}
	Rprintf("\n");
}*/
/*Rprintf("OK2\n");*/
							/* here the new partition is evaluated*/
							critFun(pM, pnr, pnc,  pnRel, pisTwoMode, pisSym,  pdiag, pnColClus, pnRowClus, ptempnUnitsRowClu, ptempnUnitsRowClu, ptemprowParArr, ptemprowParArr, papproaches, pmaxBlockTypes, pnBlockTypeByBlock, pblocks, ptempIM, ptempEM, ptempEarr, ptemperr, pjustChange, prowCluChange, prowCluChange, psameIM, pregFun, phomFun, pusePreSpec,  ppreSpecM, pcombWeights, pmulReg, pnrInSetByClusters);
/*Rprintf("Error after move = %.2f\n", *ptemperr);*/
/*Rprintf("OK3\n");*/
							if (*ptemperr< (*perr)) {
/*Rprintf("OK4a\n");*/
/*Rprintf("################################################################\n");								*/
								improve++;
								
								if (*ptemperr< (*pbesterr)) {
									updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr);
									*psameErr = 1;
								}
								if (improve >= (*pnImp)) break;
							} 
							
							if (*ptemperr == (*pbesterr)) {*psameErr += 1;}
/*Rprintf("OK4b\n");*/

							
							/* undo change found */
							ptempnUnitsRowClu[iClu2]--;	/* this line must be before the line below */
							ptemprowParArr[iClu2*(*pnr)+ptempnUnitsRowClu[iClu2]] = prowParArr[iClu2*(*pnr)+ptempnUnitsRowClu[iClu2]];
							ptemprowParArr[iClu*(*pnr)+iUnit]=prowParArr[iClu*(*pnr)+iUnit];
							ptempnUnitsRowClu[iClu]++; /* this line must be after the above line */

							/* temp values must be set to equal permanent to be updated as needed if justChange is used*/
							if(*pjustChange){
								/* temp result - copy "regular" to temp - start*/
								/* image matrix */
								for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
									ptempIM[i] = pIM[i];
								}

								/* error matrix */
								for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
									ptempEM[i] = pEM[i];
								}

								/* error array by block types*/
								for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
									ptempEarr[i] = pEarr[i];
								}
								/* temp result  - end*/
							}
						}
/*Rprintf("OK5\n");*/
						/*check the exchange of units only if iClu1 < iClu2 to avoid repeating the same move*/
						if((iClu < iClu2)&&(*pjustMove==0)){
							/* to make the order of evaluation random - start*/
							int rndUnitsInClu2[pnUnitsRowClu[iClu2]];
							for(int i=0;i<pnUnitsRowClu[iClu2];i++){
								rndUnitsInClu2[i]=i;
							}
							/* to make the order of evaluation random - end*/

							for(int iRndUnit2=0;iRndUnit2 < pnUnitsRowClu[iClu2];iRndUnit2++){
/*Rprintf("Start loop unit in cluster 2\n");*/
								/* to make the order of evaluation random - start*/
								rnd=randomInt(pnUnitsRowClu[iClu2]-iRndUnit2);
/*Rprintf("rnd = %i, pnUnitsRowClu[iClu2]-iRndUnit2= %i\n", rnd, pnUnitsRowClu[iClu2]-iRndUnit2);*/

								iUnit2=rndUnitsInClu2[rnd];
/*Rprintf("rndUnitsInClu2[rnd] = %i\n", rndUnitsInClu2[rnd]);						*/
								rndUnitsInClu2[rnd]=rndUnitsInClu2[pnUnitsRowClu[iClu2]-iRndUnit2-1];
/*Rprintf("rndUnitsInClu2[rnd] = %i\n", rndUnitsInClu2[rnd]);*/

								/* to make the order of evaluation random - end*/
								int unit2=ptemprowParArr[iClu2*(*pnr)+iUnit2];
								ptemprowParArr[iClu2*(*pnr)+iUnit2]=ptemprowParArr[iClu*(*pnr)+iUnit];
								ptemprowParArr[iClu*(*pnr)+iUnit]=unit2;

								/* here the new partition is evaluated*/
/*Rprintf("OK2-2\n");*/
/*Rprintf("iClu = %i, iClu2= %i, iUnit=%i, iUnit2=%i\n", iClu, iClu2, iUnit, iUnit2);*/
/*Rprintf("nClu = %i, nCluOld= %i, nClu2 = %i, nCluOld2= %i\n", ptempnUnitsRowClu[iClu], pnUnitsRowClu[iClu], ptempnUnitsRowClu[iClu2], pnUnitsRowClu[iClu2]);*/
/*for(int i1=0;i1<(*pnRowClus);i1++){*/
/*	Rprintf("cluster = %i, unitsCluster= %i: ", i1, ptempnUnitsRowClu[i1]);*/
/*	for(int i2=0;i2<(ptempnUnitsRowClu[i1]);i2++){*/
/*		Rprintf("%i ", ptemprowParArr[i1*(*pnr)+i2]);*/
/*	}*/
/*	Rprintf("\n");*/
/*}*/


								critFun(pM, pnr, pnc,  pnRel, pisTwoMode, pisSym,  pdiag, pnColClus, pnRowClus, ptempnUnitsRowClu, ptempnUnitsRowClu, ptemprowParArr, ptemprowParArr, papproaches, pmaxBlockTypes, pnBlockTypeByBlock, pblocks, ptempIM, ptempEM, ptempEarr, ptemperr, pjustChange, prowCluChange, prowCluChange, psameIM, pregFun, phomFun, pusePreSpec,  ppreSpecM, pcombWeights, pmulReg, pnrInSetByClusters);
/*Rprintf("OK3-2\n");*/
/* Rprintf("Error after exchange = %.2f\n", *ptemperr); */

								if (*ptemperr< (*perr)) {
	/*Rprintf("OK4a\n");*/
	/*Rprintf("################################################################\n");								*/
									improve++;
									
									if (*ptemperr< (*pbesterr)) {
										updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr);
										*psameErr = 1;
									}
									if (improve >= (*pnImp)) break;
								} 
								
								if (*ptemperr == (*pbesterr)) {*psameErr += 1;}
	/*Rprintf("OK4b\n");*/

								/* undo change found */
								ptemprowParArr[iClu*(*pnr)+iUnit]=ptemprowParArr[iClu2*(*pnr)+iUnit2];
								ptemprowParArr[iClu2*(*pnr)+iUnit2]=unit2;



								/* temp values must be set to equal permanent to be updated as needed if justChange is used*/
								if(*pjustChange){
									/* temp result - copy "regular" to temp - start*/
									/* image matrix */
									for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
										ptempIM[i] = pIM[i];
									}

									/* error matrix */
									for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
										ptempEM[i] = pEM[i];
									}

									/* error array by block types*/
									for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
										ptempEarr[i] = pEarr[i];
									}
									/* temp result  - end*/
								}
								
							}								
							/*Rprintf("OK4b-2\n");*/
						}
						/*Rprintf("OK6\n");*/
						if(improve >= (*pnImp)) break;
					}
					if(improve >= (*pnImp)) break;
				}
				if(improve >= (*pnImp)) break;
			}
		}

		free(pbestnUnitsRowClu);
		free(pbestrowParArr);
		free(pbestIM);
		free(pbestEM);
		free(pbestEarr);
		free(pbesterr);

		free(ptempnUnitsRowClu);
		free(ptemprowParArr);
		free(ptempIM);
		free(ptempEM);
		free(ptempEarr);
		free(ptemperr);

	}
	PutRNGstate(); /* Write .Random.seed in R */
/*	Rprintf("test2");*/

}




















void optParMulti(const double *pM, const int *pnr, const int *pnc,  const int *pnRel, const int *pisTwoMode, const int *pisSym, const int *pdiag, const int *pjustMove, const int *pnColClus, const int *pnRowClus, int *pnUnitsRowClu, int *pnUnitsColClu, int *prowPar, int *pcolPar, int *prowParArr, int *pcolParArr,const int *papproaches, const int *pmaxBlockTypes,const int *pnBlockTypeByBlock, const int *pblocks, int *pIM, double *pEM, double *pEarr, double *perr, const int *pjustChange, int *prowCluChange, int *pcolCluChange, const int *psameIM, const int *pregFun, const int *phomFun, const int *pusePreSpec, const double *ppreSpecM, const int *pminUnitsRowCluster, const int *pminUnitsColCluster, const int *pmaxUnitsRowCluster, const int *pmaxUnitsColCluster, int *psameErr, int *pnIter, const double *pcombWeights, const int *pexchageClusters, const int *pmaxPar, int *pbestColParMatrix, int *pbestRowParMatrix, const int *pmulReg, const int *pnrInSetByClusters){
	/*
	double *pM - pointer to array or matrix representiing the (multirelational) network
	int *pnr - pointer to the number of rows
	int *pnc - pointer to the number of columns
	int *pisTwoMode - pointer to 0 (false) or 1 (true) specifying it the network is two-mode
	int *pisSym - pointer to array of length (nRel - number of relation) specifying if the matrix (for each relation) is symmetric) 0 (false) or 1 (true) 
	int *pdiag - pointer to array of length (nRel - number of relation) specifying how to treat the diagonal elements (0 - as any other value, 1 - separately, 2 - ignore)
	int *pjustMove - pointer weather to try just moves (1) or moves and exchanges (0)	
	int *pnRel - pointer to the number of relations
	int *pnColClus - pointer to the number of column clusters
	int *pnRowClus - pointer to the number of column clusters
	int *pnUnitsRowClu - pointer to the array of the number of members of each row cluster
	int *pnUnitsColClu - pointer to the array of the number of members of each col cluster
	int *prowParArr - pointer to the array of arrays (one for each row cluster) of members of each row cluster
	int *pcolParArr - pointer to the array of arrays (one for each col cluster) of members of each col cluster
	int *papproaches - pointer to the array specifying approach - one for each relation
	int *pmaxBlockTypes - pointer to maximum number of used block types
	int *pnBlockTypeByBlock - pointer to 3d array (Rel, row, col) specifying the number of used allowed block types
	int *pblocks - pointer to the 4d array (nBlockTypesByBlock, Rel, row, col) specifying allowed block types
	int *pIM - pointer to 3d array (Rel, row, col) specifying the image matrix
	double *pEM - pointer to 3d array (Rel, row, col) specifying the error for each block
	double *pEarr - pointer to the 4d array ((*pmaxBlockTypes), Rel, row, col) specifying the errors for each allowed block type - it is very important that the value is Infinitive for block types that are not allowed
	double *perr - pointer to the total error
	int *pjustChange - pointer to a value specifying if only the errors for changed clusters should be computed
	int *prowCluChange - pointer to an array holding the two row clusters where the change occurred
	int *pcolCluChange - pointer to an array holding the col row clusters where the change occurred
	int *psameIM - pointer to 0 (false) or 1 (true) specifying if the image has to be the same for all relations
	int *pregFun - pointer to the 4d array ((*pmaxBlockTypes), Rel, row, col) specifying the "summary" function used in f-regular line blocks
	int *phomFun - pointer to the array (one value for each rel) function used used for computing measure of variability in sum of squares blockmodeling
	int *pusePreSpec - pointer to 4d array ((*pmaxBlockTypes), Rel, row, col) specifying weather a the pre-specified value should be used when computing inconsistency
	double *ppreSpecM - pointer to 4d array ((*pmaxBlockTypes), Rel, row, col) specifying the pre-specified value to be used when computing inconsistency
	int *pminUnitsRowCluster - pointer to the minimum number of units in row cluster
	int *pminUnitsColCluster - pointer to the minimum number of units in col cluster
	int *pmaxUnitsRowCluster - pointer to the maximum number of units in row cluster
	int *pmaxUnitsColCluster - pointer to the maximum number of units in col cluster
	double *pcombWeights - pointer to a array of weights of the same dimensions as blocks
	int *pexchageClusters - pointer to a matrix (nRowClust, nColClus) showing which clusters are exchangeable
	int *pmaxPar - pointer to maximum number of "best" partitions to be saved
	int *pbestColParMatrix - pointer to maximum od pmaxPar best column partitions in a matrix
	int *pbestRowParMatrix - pointer to maximum od pmaxPar best row partitions in a matrix
	int *pnrInSetByClusters - pointer to the vector of sizes of sets for row clusters. The length of the vector is equal to the number of row clusters and gives the number of units in a set to which a given row cluster belongs to. Should be set to vector of 0s if ll for group sizes should not be computed. 	
	*/

/*Rprintf("OptParC\n");*/
/**/
/*Rprintf("pM: ");*/
/*for( int i=0;i<(*pnr)*(*pnc)*(*pnRel);i++){*/
/*	Rprintf("%f ", pM[i]);*/
/*}*/
/*Rprintf("\n");*/

/*	int *pzero;
	pzero = (int *) malloc(sizeof(int));
	*pzero = 0; */
	int zero = 0;
	int rnd;

/*	Rprintf("test1");*/

	GetRNGstate();	/* Get .Random.seed from R */
	if(*pisTwoMode){
		Rprintf("Optimization of two-mode networks is corrently supported through one-mode networks.\n");
	} else {

		critFun(pM, pnr, pnc,  pnRel, pisTwoMode, pisSym,  pdiag, pnColClus, pnRowClus, pnUnitsRowClu, pnUnitsColClu, prowParArr, pcolParArr, papproaches, pmaxBlockTypes, pnBlockTypeByBlock, pblocks, pIM, pEM, pEarr, perr, &zero, prowCluChange, pcolCluChange, psameIM, pregFun, phomFun, pusePreSpec,  ppreSpecM, pcombWeights, pmulReg, pnrInSetByClusters);
/*Rprintf("Initial error = %.2f\n", *perr);*/

		/* prepare temporary objects - start*/


		/* best result  - start*/

		/* partition*/
		int *pbestrowParArr;
		int *pbestnUnitsRowClu;
		pbestnUnitsRowClu = (int *) malloc((*pnRowClus)*sizeof(int));
		pbestrowParArr = (int *) malloc((*pnRowClus)*(*pnc)*sizeof(int));
		for(int i=0;i<*pnRowClus;i++){
			pbestnUnitsRowClu[i] = pnUnitsRowClu[i];
		}
		for(int i=0;i<((*pnRowClus)*(*pnc));i++){
			pbestrowParArr[i] = prowParArr[i];
		}

		int *pbestrowPar= (int *) malloc((*pnc)*sizeof(int));

/* the next 3 lines are not necesarry */
		for(int i=0;i<(*pnc);i++){
			pbestrowPar[i] = prowPar[i];
		}

		/* image matrix */
		int *pbestIM;
		pbestIM = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestIM[i] = pIM[i];
		}

		/* number of block types by block - not needed
		int *pbestnBlockTypeByBlock;
		pbestnBlockTypeByBlock = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestnBlockTypeByBlock[i] = pnBlockTypeByBlock[i];
		} */


		/* error matrix */
		double *pbestEM;
		pbestEM = (double *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestEM[i] = pEM[i];
		}

		/* error array by block types*/
		double *pbestEarr;
		pbestEarr = (double *) malloc((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			pbestEarr[i] = pEarr[i];
		}


		double *pbesterr;
		pbesterr = (double *) malloc(sizeof(double));
		*pbesterr = *perr;

		/* best result  - end*/




		/* temp result  - start*/

		/* partition*/
		int *ptemprowParArr;
		int *ptempnUnitsRowClu;
		ptempnUnitsRowClu = (int *) malloc((*pnRowClus)*sizeof(int));
		ptemprowParArr = (int *) malloc((*pnRowClus)*(*pnc)*sizeof(int));
		for(int i=0;i<*pnRowClus;i++){
			ptempnUnitsRowClu[i] = pnUnitsRowClu[i];
		}
		for(int i=0;i<((*pnRowClus)*(*pnc));i++){
			ptemprowParArr[i] = prowParArr[i];
		}

		int *ptemprowPar= (int *) malloc((*pnc)*sizeof(int));

		/* image matrix */
		int *ptempIM;
		ptempIM = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempIM[i] = pIM[i];
		}

		/* number of block types by block - not needed
		int *ptempnBlockTypeByBlock;
		ptempnBlockTypeByBlock = (int *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(int));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempnBlockTypeByBlock[i] = pnBlockTypeByBlock[i];
		} */


		/* error matrix */
		double *ptempEM;
		ptempEM = (double *) malloc((*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempEM[i] = pEM[i];
		}

		/* error array by block types*/
		double *ptempEarr;
		ptempEarr = (double *) malloc((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus)*sizeof(double));
		for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
			ptempEarr[i] = pEarr[i];
		}

		double *ptemperr;
		ptemperr = (double *) malloc(sizeof(double));
		*ptemperr = *perr;

		/*double *ptempBestErr;
		ptempBestErr = (double *) malloc(sizeof(double));
		*ptempBestErr = *perr;*/


		/* temp result  - end*/



		/* prepare temporary objects - end*/



		int improve=1;
/*Rprintf("OK1\n");*/
		*psameErr = 1;

		/* loop until no improvement is found */
		*pnIter=0;

		while(improve){
			*pnIter = *pnIter + 1;
			/* copy best results to permanent  - start*/
			/* partition*/
			for(int i=0;i<*pnRowClus;i++){
				pnUnitsRowClu[i] = pbestnUnitsRowClu[i];
			}
			for(int i=0;i<((*pnRowClus)*(*pnc));i++){
				prowParArr[i] = pbestrowParArr[i];
			}

			/* image matrix */
			for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
				pIM[i] = pbestIM[i];
			}

			/* error matrix */
			for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
				pEM[i] = pbestEM[i];
			}

			/* error array by block types*/
			for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
				pEarr[i] = pbestEarr[i];
			}

			*perr = *pbesterr;

			/* copy best results to permanent  - end*/


			improve=0;

			/*update temp results */
			updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr);

/*Rprintf("nIter = %i\n", *pnIter);*/
/*Rprintf("pbestrowPar: ");*/
/*for( int i=0;i<*pnc;i++){*/
/*	Rprintf("%i ", pbestrowPar[i]);*/
/*}*/
/*Rprintf("\n");*/


			/* to make the order of evaluation random - start */
			/*  randomize(); does not work */
/*			int rnd;
			int rndClusters[*pnRowClus];
			for(int i=0;i<*pnRowClus;i++){
				rndClusters[i]=i;
			}	*/
			/* to make the order of evaluation random - end */

/*			int iClu, iClu2, iUnit, iUnit2;*/

			/* a loop over all clusters - random order */
			for(int iClu=0;iClu<*pnRowClus;iClu++){
/*Rprintf("Start loop cluster 1\n");*/

				/* to make the order of evaluation random - start*/
/*				rnd=randomInt(*pnRowClus-iRndClu);
				iClu=rndClusters[rnd];

				rndClusters[rnd]=rndClusters[*pnRowClus-iRndClu-1]; */
				/* to make the order of evaluation random - end*/
				prowCluChange[0]=iClu;

				/* a loop over all units inside clusters*/

				/* to make the order of evaluation random - start*/
/*				int rndUnitsInClu[pnUnitsRowClu[iClu]];
				for(int i=0;i<pnUnitsRowClu[iClu];i++){
					rndUnitsInClu[i]=i;
				}	*/
				/* to make the order of evaluation random - end*/

/*Rprintf("rndUnitsInClu: ");*/
/*for( int i=0;i<pnUnitsRowClu[iClu];i++){*/
/*	Rprintf("%i ", rndUnitsInClu[i]);*/
/*}*/
/*Rprintf("\n");*/

				for(int iUnit=0;iUnit < pnUnitsRowClu[iClu];iUnit++){
/*				for(int iRndUnit=0;iRndUnit < pnUnitsRowClu[iClu];iRndUnit++){					*/
/*Rprintf("Start loop unit in cluster 1\n");*/
					/* to make the order of evaluation random - start */
/*					rnd=randomInt(pnUnitsRowClu[iClu]-iRndUnit);
					iUnit=rndUnitsInClu[rnd];
					rndUnitsInClu[rnd]=rndUnitsInClu[pnUnitsRowClu[iClu] - iRndUnit - 1];
*/
/*Rprintf("OK 1.03\n");*/
/*Rprintf("rnd = %i, pnUnitsRowClu[iClu]-iRndUnit= %i, iUnit=%i, rndUnitsInClu[rnd]=%i\n", rnd, pnUnitsRowClu[iClu]-iRndUnit, iUnit, rndUnitsInClu[rnd]);			*/
/*Rprintf("Selected unit ID = %i\n", ptemprowParArr[iClu*(*pnr)+iUnit]);*/
/*Rprintf("rndUnitsInClu: ");*/
/*for( int i=0;i<pnUnitsRowClu[iClu];i++){*/
/*	Rprintf("%i ", rndUnitsInClu[i]);*/
/*}*/
/*Rprintf("\n");*/

					/* to make the order of evaluation random - end*/


					/* to make the order of evaluation random - start*/
/*					int rndClusters2[*pnRowClus-1];
					for(int i=0;i<*pnRowClus;i++){
						if(i < iClu) {
							rndClusters2[i]=i;
						} else if (i>iClu){
							rndClusters2[i-1]=i;
						}
					}
*/
/*Rprintf("OK 1.04\n");*/

					/* to make the order of evaluation random - end*/

					/* a loop over all other clusters - random order */
/*					for(int iRndClu2=0;iRndClu2<(*pnRowClus-1);iRndClu2++){		*/
					for(int iClu2=0;iClu2<(*pnRowClus);iClu2++){
						if(iClu==iClu2) continue;
/*Rprintf("Start loop cluster 2\n");*/
						/* to make the order of evaluation random - start*/
/*						rnd=randomInt(*pnRowClus - 1 - iRndClu2);	*/
/*Rprintf("rnd = %i, *pnRowClus - 1 - iRndClu2= %i\n", rnd, *pnRowClus - 1 - iRndClu2);			*/
/*						iClu2=rndClusters2[rnd];
						rndClusters2[rnd]=rndClusters2[*pnRowClus - 2 - iRndClu2]; 		*/
/*Rprintf("rndClusters2[rnd] = %i\n", rndClusters2[rnd]);*/
						prowCluChange[1]=iClu2;
/*Rprintf("Test exchange - start\n");*/
						if (!pexchageClusters[iClu*(*pnRowClus)+iClu2]){
							continue;
						}
/*Rprintf("Test exchange - end\n");*/

						/* to make the order of evaluation random - end*/

						if((pnUnitsRowClu[iClu]>(*pminUnitsRowCluster))&&(pnUnitsRowClu[iClu2]<(*pmaxUnitsRowCluster))){
/*Rprintf("OK1.1\n");*/
							/* move unit to another cluster */
							ptemprowParArr[iClu2*(*pnr)+ptempnUnitsRowClu[iClu2]]=ptemprowParArr[iClu*(*pnr)+iUnit];
/*Rprintf("OK1.2\n");*/

							ptempnUnitsRowClu[iClu2]++;	/* this line must be after the above line */
/*Rprintf("OK1.3\n");*/
							ptempnUnitsRowClu[iClu]--;	/* this line must be before the line below */
/*Rprintf("OK1.4\n");*/
							ptemprowParArr[iClu*(*pnr)+iUnit]=ptemprowParArr[iClu*(*pnr)+ptempnUnitsRowClu[iClu]];

/*Rprintf("iClu = %i, iClu2= %i, iUnit=%i\n", iClu, iClu2, iUnit);*/
/*Rprintf("nClu = %i, nCluOld= %i, nClu2 = %i, nCluOld2= %i\n", ptempnUnitsRowClu[iClu], pnUnitsRowClu[iClu], ptempnUnitsRowClu[iClu2], pnUnitsRowClu[iClu2]);*/
/*Rprintf("prowCluChange: %i, %i \n", prowCluChange[0], prowCluChange[1]);*/
/*for(int i1=0;i1<(*pnRowClus);i1++){*/
/*	Rprintf("cluster = %i, unitsCluster= %i: ", i1, ptempnUnitsRowClu[i1]);*/
/*	for(int i2=0;i2<(ptempnUnitsRowClu[i1]);i2++){*/
/*		Rprintf("%i ", ptemprowParArr[i1*(*pnr)+i2]);*/
/*	}*/
/*	Rprintf("\n");*/
/*}*/
/*Rprintf("OK2\n");*/
							/* here the new partition is evaluated*/
							critFun(pM, pnr, pnc,  pnRel, pisTwoMode, pisSym,  pdiag, pnColClus, pnRowClus, ptempnUnitsRowClu, ptempnUnitsRowClu, ptemprowParArr, ptemprowParArr, papproaches, pmaxBlockTypes, pnBlockTypeByBlock, pblocks, ptempIM, ptempEM, ptempEarr, ptemperr, pjustChange, prowCluChange, prowCluChange, psameIM, pregFun, phomFun, pusePreSpec,  ppreSpecM, pcombWeights, pmulReg, pnrInSetByClusters);
/*Rprintf("Error after move = %.2f\n", *ptemperr);*/
/*Rprintf("Error array and blocks:\n");*/
/*int ind2d, ind3d, ind4d;*/
/*for(int iColClu=0;iColClu<*pnColClus;iColClu++){*/
/*	Rprintf("\niColClu = %i\n", iColClu);*/
/*	for(int iRowClu=0;iRowClu<*pnRowClus;iRowClu++){*/
/*		Rprintf("iRowClu = %i\n", iRowClu);*/
/*		ind2d=iColClu*(*pnRowClus) + iRowClu;*/
/*		for(int iRel=0; iRel<(*pnRel);iRel++){*/
/*			Rprintf("iRel = %i:\n", iRel);*/
/*			ind3d= (ind2d*(*pnRel)+ iRel);*/
/*			for(int iBlockType=0;iBlockType<(pnBlockTypeByBlock[ind3d]);iBlockType++){*/
/*				ind4d=ind3d*(*pmaxBlockTypes)+iBlockType;*/
/*				Rprintf("Blocktype = %i, err = %.5f \n", pblocks[ind4d], ptempEarr[ind4d]);*/
/*			}*/
/*		}*/
/*	}*/
/*}*/


/*Rprintf("OK3\n");*/
							if (*ptemperr< (*pbesterr)) {
/*								Rprintf("Error after move = %.2f\n", *ptemperr);*/
								*psameErr=1;
								*pbesterr= *ptemperr;

								updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr);

								parArr2Vec(pnc, pnRowClus, ptempnUnitsRowClu, ptemprowParArr, pbestrowPar);
								for(int i=0;i<(*pnc);i++){
									pbestRowParMatrix[i] = pbestrowPar[i];
								}

								/* Zdajle poskusam narediti tako, da bo program sel cez vsa mozna razbitja in shranil doloceno stevilo najboljsih
								Torej da se zanka ne bo zakljucila, ko se bo naslo prvo boljse razbitje
								Pazi da bos popravil spremembe, tako tko spodaj, na zactku iteracije pa jih je potrebno ponovno udejanjiti!!!

								Mogoce se da kako bolje to narediti!!!
								*/


/*Rprintf("OK4a\n");*/
/*Rprintf("################################################################\n");*/
								improve=1;
							} else {
								if (*ptemperr == (*pbesterr)) {
									*psameErr += 1;

									int randTemp=randomInt(*psameErr);
/*									Rprintf("Error after move = %.2f\n", *ptemperr);*/
/*									Rprintf("rndUpdate = %i\n", randTemp);*/
									if(randTemp == 0){
										updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr);

										parArr2Vec(pnc, pnRowClus, ptempnUnitsRowClu, ptemprowParArr, pbestrowPar);

										if(*psameErr <= *pmaxPar){
											for(int i=0;i<(*pnc);i++){
												pbestRowParMatrix[((*psameErr)-1)*(*pnc)+i] = pbestrowPar[i];
											}
										}else{
											rnd=randomInt(*psameErr);
/*											Rprintf("rndOverwrite = %i\n", rnd);*/
											if (rnd< *pmaxPar){
												for(int i=0;i<(*pnc);i++){
													pbestRowParMatrix[rnd*(*pnc)+i] = pbestrowPar[i];
												}
											}
										}
									} else{
										parArr2Vec(pnc, pnRowClus, ptempnUnitsRowClu, ptemprowParArr, ptemprowPar);

										if(*psameErr <= *pmaxPar){
											for(int i=0;i<(*pnc);i++){
												pbestRowParMatrix[((*psameErr)-1)*(*pnc)+i] = ptemprowPar[i];
											}
										}else{
											rnd=randomInt(*psameErr);
/*											Rprintf("Error after move = %.2f\n", *ptemperr);*/
/*											Rprintf("rndOverwrite = %i\n", rnd);*/
											if (rnd< *pmaxPar){
												for(int i=0;i<(*pnc);i++){
													pbestRowParMatrix[rnd*(*pnc)+i] = ptemprowPar[i];
												}
											}
										}
									}
								}
							}
/*Rprintf("OK4b\n");*/

							/* undo change found */
							ptempnUnitsRowClu[iClu2]--;	/* this line must be before the line below */
							ptemprowParArr[iClu2*(*pnr)+ptempnUnitsRowClu[iClu2]] = prowParArr[iClu2*(*pnr)+ptempnUnitsRowClu[iClu2]];
							ptemprowParArr[iClu*(*pnr)+iUnit]=prowParArr[iClu*(*pnr)+iUnit];
							ptempnUnitsRowClu[iClu]++; /* this line must be after the above line */

							/* temp values must be set to equal permanent to be updated as needed if justChange is used*/
							if(*pjustChange){
								/* temp result - copy "regular" to temp - start*/
								/* image matrix */
								for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
									ptempIM[i] = pIM[i];
								}

								/* error matrix */
								for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
									ptempEM[i] = pEM[i];
								}

								/* error array by block types*/
								for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
									ptempEarr[i] = pEarr[i];
								}
								/* temp result  - end*/
							}

/*Rprintf("OK5\n");*/
						}

						/*check the exchange of units only if iClu1 < iClu2 to avoid repeating the same move */
						if((iClu < iClu2)&&(*pjustMove==0)){
							/* to make the order of evaluation random - start*/
/*							int rndUnitsInClu2[pnUnitsRowClu[iClu2]];
							for(int i=0;i<pnUnitsRowClu[iClu2];i++){
								rndUnitsInClu2[i]=i;
							}		*/
							/* to make the order of evaluation random - end*/

/*							for(int iRndUnit2=0;iRndUnit2 < pnUnitsRowClu[iClu2];iRndUnit2++){		*/
							for(int iUnit2=0;iUnit2 < pnUnitsRowClu[iClu2];iUnit2++){
/*Rprintf("Start loop unit in cluster 2\n");*/
								/* to make the order of evaluation random - start*/
/*								rnd=randomIntpnUnitsRowClu[iClu2]-iRndUnit2);		*/
/*Rprintf("rnd = %i, pnUnitsRowClu[iClu2]-iRndUnit2= %i\n", rnd, pnUnitsRowClu[iClu2]-iRndUnit2);*/

/*								iUnit2=rndUnitsInClu2[rnd];
								rndUnitsInClu2[rnd]=rndUnitsInClu2[pnUnitsRowClu[iClu2]-iRndUnit2-1];		*/
/*Rprintf("rndUnitsInClu2[rnd] = %i\n", rndUnitsInClu2[rnd]);						*/

								/* to make the order of evaluation random - end*/
								int unit2=ptemprowParArr[iClu2*(*pnr)+iUnit2];
								ptemprowParArr[iClu2*(*pnr)+iUnit2]=ptemprowParArr[iClu*(*pnr)+iUnit];
								ptemprowParArr[iClu*(*pnr)+iUnit]=unit2;

								/* here the new partition is evaluated*/
/*Rprintf("OK2-2\n");*/
/*Rprintf("iClu = %i, iClu2= %i, iUnit=%i, iUnit2=%i\n", iClu, iClu2, iUnit, iUnit2);*/
/*Rprintf("nClu = %i, nCluOld= %i, nClu2 = %i, nCluOld2= %i\n", ptempnUnitsRowClu[iClu], pnUnitsRowClu[iClu], ptempnUnitsRowClu[iClu2], pnUnitsRowClu[iClu2]);*/
/*Rprintf("prowCluChange: %i, %i \n", prowCluChange[0], prowCluChange[1]);*/
/*for(int i1=0;i1<(*pnRowClus);i1++){*/
/*	Rprintf("cluster = %i, unitsCluster= %i: ", i1, ptempnUnitsRowClu[i1]);*/
/*	for(int i2=0;i2<(ptempnUnitsRowClu[i1]);i2++){*/
/*		Rprintf("%i ", ptemprowParArr[i1*(*pnr)+i2]);*/
/*	}*/
/*	Rprintf("\n");*/
/*}*/


								critFun(pM, pnr, pnc,  pnRel, pisTwoMode, pisSym,  pdiag, pnColClus, pnRowClus, ptempnUnitsRowClu, ptempnUnitsRowClu, ptemprowParArr, ptemprowParArr, papproaches, pmaxBlockTypes, pnBlockTypeByBlock, pblocks, ptempIM, ptempEM, ptempEarr, ptemperr, pjustChange, prowCluChange, prowCluChange, psameIM, pregFun, phomFun, pusePreSpec,  ppreSpecM, pcombWeights, pmulReg, pnrInSetByClusters);
/*Rprintf("OK3-2\n");*/
/* Rprintf("Error after exchange = %.2f\n", *ptemperr);*/

								if (*ptemperr< (*pbesterr)) {
/*									Rprintf("Error after exchange = %.2f\n", *ptemperr);*/
									*psameErr=1;
									*pbesterr= *ptemperr;

									updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr);

									parArr2Vec(pnc, pnRowClus, ptempnUnitsRowClu, ptemprowParArr, pbestrowPar);
									for(int i=0;i<(*pnc);i++){
										pbestRowParMatrix[i] = pbestrowPar[i];
									}

									improve=1;
								} else {
									if (*ptemperr == (*pbesterr)) {
										*psameErr += 1;

										int randTemp=randomInt(*psameErr);
/*										Rprintf("Error after exchange = %.2f\n", *ptemperr);*/
/*										Rprintf("rndUpdate = %i\n", randTemp);*/
										if(randTemp == 0){
											updateResults(pnc, pnRel, pnColClus, pnRowClus, pmaxBlockTypes, ptempnUnitsRowClu, ptemprowParArr, ptempIM, ptempEM, ptempEarr, ptemperr, pbestnUnitsRowClu, pbestrowParArr, pbestIM, pbestEM, pbestEarr, pbesterr);

											parArr2Vec(pnc, pnRowClus, ptempnUnitsRowClu, ptemprowParArr, pbestrowPar);

											if(*psameErr <= *pmaxPar){
												for(int i=0;i<(*pnc);i++){
													pbestRowParMatrix[((*psameErr)-1)*(*pnc)+i] = pbestrowPar[i];
												}
											}else{
												rnd=randomInt(*psameErr);
/*												Rprintf("Error after exchange = %.2f\n", *ptemperr);*/
/*												Rprintf("rndOverwrite = %i\n", rnd);*/
												if (rnd< *pmaxPar){
													for(int i=0;i<(*pnc);i++){
														pbestRowParMatrix[rnd*(*pnc)+i] = pbestrowPar[i];
													}
												}
											}
										} else{
											parArr2Vec(pnc, pnRowClus, ptempnUnitsRowClu, ptemprowParArr, ptemprowPar);

											if(*psameErr <= *pmaxPar){
												for(int i=0;i<(*pnc);i++){
													pbestRowParMatrix[((*psameErr)-1)*(*pnc)+i] = ptemprowPar[i];
												}
											}else{
												rnd=randomInt(*psameErr);
/*												Rprintf("Error after exchange = %.2f\n", *ptemperr);*/
/*												Rprintf("rndOverwrite = %i\n", rnd);*/
												if (rnd< *pmaxPar){
													for(int i=0;i<(*pnc);i++){
														pbestRowParMatrix[rnd*(*pnc)+i] = ptemprowPar[i];
													}
												}
											}
										}

									}
								}
								ptemprowParArr[iClu*(*pnr)+iUnit]=ptemprowParArr[iClu2*(*pnr)+iUnit2];
								ptemprowParArr[iClu2*(*pnr)+iUnit2]=unit2;

								/* temp values must be set to equal permanent to be updated as needed if justChange is used*/
								if(*pjustChange){
									/* temp result - copy "regular" to temp - start*/
									/* image matrix */
									for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
										ptempIM[i] = pIM[i];
									}

									/* error matrix */
									for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
										ptempEM[i] = pEM[i];
									}

									/* error array by block types*/
									for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
										ptempEarr[i] = pEarr[i];
									}
									/* temp result  - end*/
/*Rprintf("OK4b-2\n");*/
								}
							}
						}
/*Rprintf("OK6\n");*/

/*						if(improve) break;*/
					}
/*					if(improve) break; */
				}
/*				if(improve) break; */
			}
/*Rprintf("Iteration %i completed, improve = %i\n", *pnIter, improve); */
		}
		for(int i=0;i<(*pnc);i++){
			prowPar[i] = pbestrowPar[i];
		}

		free(pbestnUnitsRowClu);
		free(pbestrowParArr);
		free(pbestIM);
		free(pbestEM);
		free(pbestEarr);
		free(pbesterr);
		free(pbestrowPar);

		free(ptempnUnitsRowClu);
		free(ptemprowParArr);
		free(ptempIM);
		free(ptempEM);
		free(ptempEarr);
		free(ptemperr);
		free(ptemprowPar);
	}
     PutRNGstate(); /* Write .Random.seed in R */
     /*
Rprintf("nIter = %i\n", *pnIter);
Rprintf("prowPar: ");
for( int i=0;i<*pnc;i++){
	Rprintf("%i ", prowPar[i]);
}
Rprintf("\n");

Rprintf("pnUnitsRowClu: ");
for(int i=0;i<*pnRowClus;i++){
	Rprintf("%i ", pnUnitsRowClu[i]);
}
Rprintf("\n");

Rprintf("prowParArr: ");
for(int i=0;i<((*pnRowClus)*(*pnc));i++){
	Rprintf("%i ", prowParArr[i]);
}
Rprintf("\n");


Rprintf("pIM: ");
for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
	Rprintf("%i ", pIM[i]);
}
Rprintf("\n");

Rprintf("pEM: ");
for(int i=0;i<((*pnRel)*(*pnRowClus)*(*pnColClus));i++){
	Rprintf("%lf ", pEM[i]);
}
Rprintf("\n");


Rprintf("pEarr: ");
for(int i=0;i<((*pmaxBlockTypes)*(*pnRel)*(*pnRowClus)*(*pnColClus));i++){
	Rprintf("%lf ", pEarr[i]);
}
Rprintf("\n");


Rprintf("perr: %lf \n", *perr);

Rprintf("prowCluChange: %i, %i \n", prowCluChange[0], prowCluChange[1]);

Rprintf("psameErr: %i \n", *psameErr);

Rprintf("pnIter: %i \n", *pnIter);

Rprintf("pbestRowParMatrix: ");
for(int i=0;i<((*pnc)*(*pmaxPar));i++){
	Rprintf("%i ", pbestRowParMatrix[i]);
}
Rprintf("\n");


Rprintf("Function completed\n");
*/

}


