/*! \file findinits.c

by Andrew Redd
find initial values for the single algorithm
*/

#include "pfda.h"
#include <R_ext/Lapack.h>



void pfdafindinits( double const * const y, double const * const z, int const * const nobs,int const * const M, int const * const N, int const * const ka,int const * const kb,double const * const B, int const * const p, double const * const delta, //input values
double* const tm,double* const tn,double* const tf, double* const tg, double* const Da, double* const Db, double* const L, double* const seps, double* const sxi, double* const Alpha, double* const Beta, double* const Saa, double* const Sab, double* const Sbb) //return values
{
/*
inputs:
y		M×1 vector of y values
z		M×1 vector of zvalues
obs		M×1 vector of subject identifiers
nobs		N×1 vector of number of observations per subject
M		number of observations
N		number of subjects
ka		number of principle components used for modeling y
kb		number of principle components used for modeling z
B		M×p matrix of timepoints evaluated at the orthogonalized basis functions
p		number of orthogonal basis functions
delta		for ridge regression type computations

outputs:	name|Description(name in R)
tm		p×1 vector estimating mean structore for y (theta_mu)
tn		p×1 vector estimating mean structore for z (theta_nu)
tf		p×ka matrix representing the subject specific portion of the estimate for y(Theta_f)
tg		p×kb matrix representing the subject specific portion of the estimate for z (Theta_g)
Da		ka vector representing the variances of Alpha (D_Alpha)
Db		kb vector representing the variance of Beta (D_Beta)
L		kb×ka matrix representing the covariance of A and B (Lambda)
seps		variance of y residuals (sigma_epsilon)
sxi		variacne of z residuals (sigma_xi)
Alpha		N×ka matrix of PC coefficients of y (Alpha)
Beta		N×kb matrix of PC coefficients of z(Beta)
Saa		ka×ka×N 3 dimentional array of variance components (Sigma_alpha_alpha)
Sab		ka×kb×N 3 dimentional array of variance components (Sigma_alpha_beta)
Sbb		kb×kb×N 3 dimentional array of variance components (Sigma_beta_beta)
*/

//find tm tn
int i,j,k; //reserved indexes for itteration.
double* BtB; //pXp matrix of the inner product of B
int svr; //results for solve opperations
BtB = Idiag(*p);// initialize BtB to identity matrix
dsyrk_(&Upper, &Trans, p, M, &dOne, B, M, delta, BtB, p);//BtB=B'B+deltaI 
dgemv_(&Trans,M,p,&dOne,B,M,y,&one,&dzero,tm,&one);  //  compute tm = B^t y
dposv_(&Upper, p, &one, BtB, p, tm, p, &svr);
if(svr>0)Rf_error("PFDA ERR1: Leading minor of order %i is not positive definite",svr); 
dgemv_(&Trans,M, p,&dOne,B,M,z,&one,&dzero,tn,&one);  //  compute tn = B^t z
dtrsv_(&Upper, &Trans, &NonUnit, p, BtB, p, tn, &one);
dtrsv_(&Upper, &NoTrans, &NonUnit, p, BtB, p, tn, &one);
//end tm tn //begin tf tg Da Db ///////////////////////////////////////////////////////////////////////////////////////
double *yc, *zc, *btbi, *stfa, *stgb;
double const *Bi, *yci, *zci;
double **tfa, **tgb;  //two dimentional makes most sense for these.
yc = (double*)R_alloc(*M,sizeof(double));for(i=0;i<*M;i++){yc[i]=y[i];}//copy y into yc
zc = (double*)R_alloc(*M,sizeof(double));for(i=0;i<*M;i++){zc[i]=z[i];}//copy z into zc
dgemv_( &NoTrans, M, p, &mOne, B, M, tm, &one, &dOne,yc,&one);// yc = yc-B*tm
dgemv_( &NoTrans, M, p, &mOne, B, M, tn, &one, &dOne,zc,&one);// zc = zc-B*tn
stfa = (double*)S_alloc((*p)*(*p+1)/2,sizeof(double));
stgb = (double*)S_alloc((*p)*(*p+1)/2,sizeof(double));  //S_alloc to initialize to zero.  packed.
Bi=B;yci=yc;zci=zc;
tfa = (double**)S_alloc((*N),sizeof(double*));
tgb = (double**)S_alloc((*N),sizeof(double*));
btbi= (double*)S_alloc((*p)*(*p),sizeof(double));//btbi is reused
for(k=0;k<*N;k++){
	//allocate arrays
	tfa[k] = (double*)S_alloc((*p),sizeof(double));
	tgb[k] = (double*)S_alloc((*p),sizeof(double));

	resetI(btbi, *p);
	
	dsyrk_( &Upper, &Trans, p, &(nobs[k]), &dOne, Bi, M, delta, btbi, p); //btbi=Bi'Bi+deltaI
	dgemv_( &Trans,&(nobs[k]),p,&dOne, Bi, M, yci,&one,&dzero,tfa[k],&one);  //  compute tfa = Bi^t yci = Bi^t (yi-Bi*tm)
	dposv_( &Upper, p, &one, btbi, p, tfa[k], p, &svr);  //(btbi)^-1 (tfa)  // Cholesky given in btbi
	if(svr>0)Rf_error("PFDA ERR2: Leading minor of order %i is not positive definite. Iteration: %d",svr,k);
	dgemv_( &Trans,&(nobs[k]),p,&dOne, Bi, M, zci,&one,&dzero,tgb[k],&one);  //  compute tgb = Bi^t zci = B^t (zi-Bi*tn)
	dtrsv_( &Upper, &Trans, &NonUnit, p, btbi, p, tgb[k],&one); // see below
	dtrsv_( &Upper, &NoTrans, &NonUnit, p, btbi, p, tgb[k], &one);  // solves  linear equation for zsince we have the cholesky.

	dspr_(  &Upper, p, &dOne, tfa[k],&one, stfa);
	dspr_(  &Upper, p, &dOne, tgb[k],&one, stgb);
	
	Bi+=nobs[k];yci+=nobs[k];zci+=nobs[k];
}
//printspmat(stgb,*p);
double *evalf, *evecf, *evalg, *evecg, *work;
evalf = (double*)R_alloc(*p,        sizeof(double));
evecf = (double*)R_alloc((*p)*(*p), sizeof(double));
evalg = (double*)R_alloc(*p,        sizeof(double));
evecg = (double*)R_alloc((*p)*(*p), sizeof(double));
work =  (double*)R_alloc(3*(*p),   	sizeof(double));
//workg = (double*)R_alloc(3*(*p),    sizeof(double));
//Compute eigen vector and eigen values
int eigenresult;
eigenresult = local_dspev(p, stfa, evalf, evecf, p, work);
if(eigenresult)Rf_error("PFDA ERR3: Could not find eigen values associated with Alpha");
eigenresult = local_dspev(p, stgb, evalg, evecg, p, work);
if(eigenresult)Rf_error("PFDA ERR4: Could not find eigen values associated with Beta");
//copy relevent vlaues to relevent places
for(i=0;i<*ka;i++){
  Da[i]=evalf[*p-i-1];
  for(j=0;j<*p;j++){
    tf[j+(*p)*i]=evecf[(*p-i-1)*(*p)+j];
  }
}
for(i=0;i<*kb;i++){
  Db[i]=evalg[*p-i-1];
  for(j=0;j<*p;j++){
    tg[j+(*p)*i]=evecg[(*p-i-1)*(*p)+j];
  }
}
//find Alpha and Beta
positivefirstrow(tf,*p,*ka);
positivefirstrow(tg,*p,*kb);

// double *ttt,*gtg;
// ttt=(double*)S_alloc(*ka*(*ka),sizeof(double));
// gtg=(double*)S_alloc(*kb*(*kb),sizeof(double));

for(k=0;k<*N;k++){
	dgemv_( &Trans,p,ka,&dOne, tf, p, tfa[k],&one,&dzero,&(Alpha[k]), N);  //  compute Alpha = tf^t tfa[i]
	dgemv_( &Trans,p,kb,&dOne, tg, p, tgb[k],&one,&dzero,&(Beta[k]), N);  //  compute Beta = tg^t tgb[i]
}
//Find Lambda
double *ata;
ata = (double*)R_alloc(*ka**ka,sizeof(double));
dsyrk_( &Upper, &Trans, ka, N, &dOne, Alpha, N, &dzero, ata, ka); //btbi=Bi'Bi+deltaI
dgemm_( &Trans, &NoTrans, ka, kb, N, &dOne, Alpha, N, Beta, N, &dzero, L, ka); //L=Alpha^t Beta
dposv_( &Upper, ka, kb, ata, ka, L, ka, &svr);  //L= (Alpha^t Alpha)^-1 Alpha^t Beta
if(svr)Rf_error("PFDA ERR5: Leading minor of order %i is not positive definite.",svr,k);

//zero out the Sigma matrices
for(i=0;i<*N**ka**ka;i++)Saa[i]=0;
for(i=0;i<*N**ka**kb;i++)Sab[i]=0;
for(i=0;i<*N**kb**kb;i++)Sbb[i]=0;

/* //Find Saa Sab and Sbb
double *suaa, *suab, *subb;
double *BTf, *BTg;
double *BTfi, *BTgi;
Bi=0;
BTfi = BTf = (double*)R_alloc(*M**ka, sizeof(double));
BTgi = BTg = (double*)R_alloc(*M**kb, sizeof(double));
dgemm_(&NoTrans, &NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, BTf, M);//compute BTf = B*Tf
dgemm_(&NoTrans, &NoTrans, M, kb, p, &dOne, B, M, tg, p, &dzero, BTg, M);//compute BTg = B*Tg
 */


}

