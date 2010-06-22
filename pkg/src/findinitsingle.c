/*! \file findinitsingle.c

by Andrew Redd
find initial values for the single algorithm
*/
#include "pfda.h"
#include <R_ext/Lapack.h>

/*! new function for computing the initial values

@MEMORY
	- dp = 2*p^2 + M + M*k+ p*N + 8*p
	- ip = 6*p
*/
void pfda_s_i(
	double       * const tm,
	double       * const tf,
	double       * const alpha,
	double       * const Da,
	double       * const aa,
	double       * const sigma,
	double const * const y, 
	int    const * const nobs,
	int    const * const M, 
	int    const * const N, 
	int    const * const k, 
	double const * const B, 
	double const * const btb,
	int    const * const p, 
	double const * const minV,
	int const * const dl,
	double * dp, int * ip)
{
	pfda_debug_cdl(debugnum_s_i){pfda_debug_msg("entering pfda_s_i\n");fflush(stdout);}
	// tm 
	double * left = pfdaAlloc_d(*p**p,&dp);
	for(int j=0;j<*p;j++)    left[j+j**p] = *minV;
	dsyrk_(&Upper, &Trans, p, M, &dOne, B, M, &dOne, left, p);
	dgemv_(&Trans,M,p,&dOne,B,M,y,&one,&dzero,tm,&one);
	pfda_debug_cdl(debugnum_s_i){
		pfda_debug_msg("tm left:\n");printmat(left,*p,*p);
		pfda_debug_msg("tm right:\n");printmat(tm,*p,one);
		fflush(stdout);
	}
	pfda_debug_cdl(debugnum_s_1_break_tm)return;
	int svr=0;
	dposv_(&Upper, p, &one, left, p, tm, p, &svr);

	// tf
	double * Ry = pfdaAlloc_d(*M,&dp);
	pfda_computeResid(Ry,y,nobs,M,N,k,B,p,tm,NULL,NULL,dl,dp);
	double * tfa = pfdaAlloc_d(*N**p,&dp);
	double * Ryi=Ry;
	double const * Bi=B;
	for(int i=0;i<*N;i++){
		for(int j=0;j<*p**p;j++) left[j] = 0;
		for(int j=0;j<*p;j++)    left[j+j**p] = *minV;
		for(int j=0;j<*p**p;j++) left[j] += btb[j+*p**p*i];
		dgemv_(&Trans,nobs+i,p,&dOne, Bi, M, Ryi,&one,&dzero,tfa+i**p,&one);
		int svr=0;
		pfda_debug_cdl(debugnum_s_i){
			pfda_debug_msg("tfa[%d] left:\n",i);printmat(left,*p,*p);
			pfda_debug_msg("tfa[%d] right:\n",i);printmat(tfa+i**p,*p,one);
			fflush(stdout);
		}
		dposv_( &Upper, p, &one, left, p, tfa+i**p, p, &svr);  //(btbi)^-1 (tfa)  // Cholesky given in btbi
		if(svr>0)pfda_error("PFDA ERR1: Leading minor of order %i is not positive definite",svr);
		else pfda_debug_cdl(debugnum_s_i){
			pfda_debug_msg("tf[%d] solution:\n",i);printmat(tfa+i**p,*p,one);
			fflush(stdout);
		}
		
		Ryi += nobs[i];
		Bi  += nobs[i];
	}
	pfda_debug_cdl(debugnum_s_i){ pfda_debug_msg("tfa:\n");printmat(tfa,*p,*N);fflush(stdout);}
	double * T=pfdaAlloc_d(*p**p,&dp);
	dsyrk_(&Upper, &NoTrans, p, N, &dOne, tfa, p, &dzero, T, p);

	pfda_debug_cdl(debugnum_s_i){
		pfda_debug_msg("tf T:\n");printmat(T,*p,*p);
		fflush(stdout);
	}
	pfda_eigens(T,	p, tf, Da, k, dl, dp, ip);
	positivefirstrow(tf,*p,*k);
	pfda_debug_cdl(debugnum_s_i_break_alpha)return;
	for(int i=0;i<*N;i++){
		dgemv_(&Trans, p,k,&dOne, tf, p, tfa+i**p, &one, &dzero, alpha+i, N);
		pfda_debug_cdl(debugnum_s_i){
			pfda_debug_msg("alpha[%d]:\n",i);printmat2(one, *k, alpha+i, N);
			fflush(stdout);
		}
		if(aa)dsyr_(&Upper, k, &dOne, alpha+i, N, aa+*k**k*i,k); 
	}
	
	if(sigma){
		pfda_computeResid(Ry,y,nobs,M,N,k,B,p,tm,tf,alpha,dl,dp);
		*sigma = ddot_(M,Ry,&one,Ry,&one);
		*sigma /= (double)*M;
	}
	pfda_debug_cdl(debugnum_s_i){pfda_debug_msg("exiting pfda_s_i\n");fflush(stdout);}
}

void pfdaFindInitSingle(
	double const * const t, 
	double const * const y, 
	int const * const nobs,
	int const * const M, 
	int const * const N, 
	int const * const ka, 
	double const * const B, 
	int const * const p, 
	double const * const delta, //input values
	double * const tm, 
	double * const tf, 
	double * const Da, 
	double * const seps, 
	double * const Alpha,
	const int * const dl) //return values
{
/*
inputs:
t		M×1 vector of time values
y		M×1 vector of y values
nobs		N×1 vector of number of observations per subject
M		number of observations
N		number of subjects
ka		number of principle components used for modeling y
B		M×p matrix of timepoints evaluated at the orthogonalized basis functions
p		number of orthogonal basis functions
delta		for ridge regression type computations

outputs:	name|Description(name in R)
tm		p×1 vector estimating mean structore for y (theta_mu)
tf		p×ka matrix representing the subject specific portion of the estimate for y(Theta_f)
Da		ka vector representing the variances of Alpha (D_Alpha)
seps		variance of y residuals (sigma_epsilon)
Alpha		N×ka matrix of PC coefficients of y (Alpha)
*/

/* Declarations*/
int i,k; //reserved indexes for itteration.
int eigenresult;
int svr; //results for solve opperations
int debug=*dl;
double *BtB; //pXp matrix of the inner product of B
double *yc=NULL, *btbi=NULL, *stfa=NULL;
const double * Bi=NULL, *yci=NULL;
double **tfa;  //two dimentional makes most sense for these.
double *eigen_values, *eigen_vectors, *work;

/* Allocation ot Memory */
yc = (double*)R_alloc(*M,sizeof(double));for(i=0;i<*M;i++){yc[i]=y[i];}//copy y into yc
stfa = (double*)S_alloc((*p)*(*p+1)/2,sizeof(double));
tfa = (double**)S_alloc((*N),sizeof(double*));
btbi= (double*)S_alloc((*p)*(*p),sizeof(double));//btbi is reused
for(k=0;k<*N;k++) tfa[k] = (double*)S_alloc((*p),sizeof(double));
eigen_values = (double*)R_alloc(*p,        sizeof(double));
eigen_vectors = (double*)R_alloc((*p)*(*p), sizeof(double));
work =  (double*)R_alloc(3*(*p),   	sizeof(double));


/* begin */
if(debug>0){Rprintf("Entering pfdaFindInitSingle\n");fflush(stdout);}
if(debug==4){
	Rprintf("Input Values:\n");
	Rprintf("M:%d\tN:%d\tka:%d\tp:%d\n",*M,*N,*ka,*p);
	Rprintf("delta: %g\n",*delta);
	Rprintf("Vector arguments:\n");
	Rprintf("t:\n");printmat(t,1,*M);
	Rprintf("y:\n");printmat(y,1,*M);
	Rprintf("nobs:\n");printmati(nobs,1,*N);
	fflush(stdout);
}
if(debug==4){
	Rprintf("Matrix Arguments:\n");
	Rprintf("B:\n");printmat(B,*M,*p);
	fflush(stdout);
}

//find tm 
if(debug>0){Rprintf("Find tm\n");fflush(stdout);}
BtB = Idiag(*p);// initialize BtB to identity matrix
dsyrk_(&Upper, &Trans, p, M, &dOne, B, M, delta, BtB, p);//BtB=B'B+deltaI 
dgemv_(&Trans,M,p,&dOne,B,M,y,&one,&dzero,tm,&one);  //  compute tm = B^t y
dposv_(&Upper, p, &one, BtB, p, tm, p, &svr);  // tm = (BtB)^-1 B^t y
if(svr>0)Rf_error("PFDA ERR1: Leading minor of order %i is not positive definite",svr);

//find tf Da 
if(debug>0){Rprintf("Find tf Da\n");fflush(stdout);}
dgemv_( &NoTrans, M, p, &mOne, B, M, tm, &one, &dOne,yc,&one);// yc = yc-B*tm
Bi=B;yci=yc;
for(k=0;k<*N;k++){
	if(debug==4){Rprintf("loop #%d\n",k);fflush(stdout);}
	resetI(btbi, *p);
	
	dsyrk_( &Upper, &Trans, p, &(nobs[k]), &dOne, Bi, M, delta, btbi, p); //btbi=Bi'Bi+deltaI
	dgemv_( &Trans,&(nobs[k]),p,&dOne, Bi, M, yci,&one,&dzero,tfa[k],&one);  //  compute tfa = Bi^t yci = Bi^t (yi-Bi*tm)
	dposv_( &Upper, p, &one, btbi, p, tfa[k], p, &svr);  //(btbi)^-1 (tfa)  // Cholesky given in btbi
	if(svr>0)Rf_error("PFDA ERR2: Leading minor of order %i is not positive definite. Iteration: %d",svr,k);

	dspr_( &Upper, p, &dOne, tfa[k],&one, stfa);
	
	Bi+=nobs[k];yci+=nobs[k];
}
//Compute eigen vector and eigen values
if(debug>0){Rprintf("Find eigenvectors and values\n");fflush(stdout);}
eigenresult = local_dspev(p, stfa, eigen_values, eigen_vectors, p, work);
if(eigenresult)Rf_error("PFDA ERR3: Could not find eigen values associated with Alpha");
//copy relevent values to relevent places
if(debug==4){
	Rprintf("eigen Vectors:\n");printmat(eigen_vectors,*p,*p);fflush(stdout);
}

int vectornum=0, cellnum=0;
for(vectornum=0;vectornum<*ka;vectornum++){
  Da[vectornum]=eigen_values[*p-vectornum-1];
  for(cellnum=0;cellnum<*p;cellnum++){
    tf[cellnum+(*p)*vectornum]=eigen_vectors[(*p-vectornum-1)*(*p)+cellnum];
  }
}
//find Alpha
if(debug>0){Rprintf("Find Alpha\n");fflush(stdout);}
positivefirstrow(tf,*p,*ka);
for(k=0;k<*N;k++){
	dgemv_( &Trans,p,ka,&dOne, tf, p, tfa[k],&one,&dzero,&(Alpha[k]), N);  //  compute Alpha = tf^t tfa[i]
}
if(debug>0){Rprintf("Exiting pfdaFindInitSingle\n");fflush(stdout);}
}
