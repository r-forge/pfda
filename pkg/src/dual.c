/*! \file dual.c

by Andrew Redd
Part of paired functional data analysis through principle component analysis.
This file performs the em algorithm for a paired response curves.
*/
#include "pfda_dual.h"
#include "pfda_error.h"
#include "pfda_gen.h"
#include "pfda_memory.h"
#include "pfda_utils.h"
#include <R_ext/Lapack.h>

void pfdaDual_m4(
	int    const * const N,
	double       * const Lambda,
	double const * const alpha,
	int    const * const ka,
	double const * const Beta,
	int    const * const kb,
	double const * const Sigma_aa,
	double const * const Sigma_ab,
	int const * const dl, double * dp)
/* DESCRIPTION
INPUTS:
	N - number of subjects
	Lambda  (kb×ka) - Matrix of regression coefficients for regression of beta on alpha
	alpha - matrix of principle component scores for Y.
	ka  - number of columns of alpha, ie. the number of principle components for Y.
	Beta - matrix of principle components scores for Z.
	kb - numbero fo columns of Beta, ie. the number of principle components for Z.
	dl - integer array for printing debugging information.
WORKING SPACE:
	work_aa - working space of at least ka^2 length
Passing On:
	Sigma_aa - variance of alpha_i
	Sigma_ab - covariance of alpha_i beta_i
OUTPUT:
	Lambda
	\[ %tex description
		\Lambda = \left( \sum\limits_{i=1}^N \beta_i\alpha_i\trans \right) \left(\sum\limits_{i=1}^N\alpha_i\alpha_i\trans\right)
	\]
MEMORY:
	dp	length =  ka^2 + ka*kb
dl:	240
*/
{ //code
	pfda_debug_cdl(debugnum_dual_m4){
		pfda_debug_msg("pfdaDual_m4 - \n");
		pfda_debug_msg("Lambda(in):\n"); printmat(Lambda, *kb, *ka);
		pfda_debug_msg("N:\n%d\n\n",*N );
		pfda_debug_msg("alpha:\n"); printmat(alpha, *N, *ka);
		pfda_debug_msg("ka:\n%d\n\n",*ka);
		pfda_debug_msg("Beta:\n");  printmat(Beta, *N, *kb);
		pfda_debug_msg("kb:\n%d\n\n",*kb);
	}
	//* pointers declarations */
	int cellnum, solve_result;
	for(cellnum=0;cellnum<(*ka**kb);cellnum++)Lambda[cellnum]=0;
	double * work_aa = pfdaAlloc_d(*ka**ka, &dp);

	pfda_sum_cond_dd(work_aa, N, alpha, Sigma_aa, ka, dl);			// work_aa += \widehat{\alpha\alpha\trans}
	pfda_sum_cond_ab(Lambda, N, alpha, ka, Beta, kb, Sigma_ab, dl);	// Lambda += \widehat{\alpha\beta\trans}

	dposv_(&Upper, ka, kb, work_aa, ka, Lambda, ka, &solve_result);
	if(solve_result!=0)pfda_error("PFDA ERR:pfdaDual_m4: Leading minor of order %i is not positive definite(Lambda)",solve_result);
	pfda_transpose(Lambda, *ka, *kb, dl, dp);
	pfdaFree(work_aa);
	pfda_debug_cdl(debugnum_dual_m4){
		pfda_debug_msg("Lambda(out):\n"); printmat(Lambda, *kb, *ka);
		pfda_debug_msg("leaving pfdaDual_m4\n");
	}
}

void pfdaDual_m5_2(
	double       * const Lambda,
	double       * const trans_f, //will hold the chol decomp after run.
	int    const * const ka,
	double const * const trans_g,
	int    const * const kb,

	int const * const dl, double*dp,int*ip)
/* DESCRIPTION
INPUTS:
	N - number of subjects
	delta - matrix of principle component scores for Y/Z.
	k  - number of columns of delta, ie. the number of principle components for Y/Z.
	Sigma_dd - variances of delta_i.
	Theta - Coefficient matrix for the principle component curves for Y/Z.
	p - the number of rows of Theta.
	dl - integer array for printing debugging information.
WORKING SPACE:
	work_dd - k^2 double
	working space for finding principle components
OUTPUT:
	delta - updated with transformation
	theta - updated with new vectors
	Transformation (k x k) - Q^TTheta
	\[ %tex description
		\Lambda <- Trans_g \Lambda Trans_f^{-1}
	\]
MEMORY:
	dp	length = ka^2 + kb*ka
	ip	length = ka
*/
{//code
	pfda_debug_cdl(debugnum_dual_m5_2){
		pfda_debug_msg("pfdaDual_m5_2 - \n");
		pfda_debug_msg("Lambda:\n"); printmat(Lambda, *kb, *ka);
		pfda_debug_msg("trans_f:\n"); printmat(trans_f, *ka, *ka);
		pfda_debug_msg("ka:\n%d\n\n",*ka);
		pfda_debug_msg("trans_g:\n"); printmat(trans_g, *kb, *kb);
		pfda_debug_msg("kb:\n%d\n\n",*kb);
		fflush(stdout);
	}
	double * inv_trans_f = pfdaAlloc_d(*ka**ka, &dp);
	int solve_info=0;
	resetI(inv_trans_f,*ka);
	int * ipiv =pfdaAlloc_i (*ka, &ip);
	dgesv_(ka,ka,trans_f, ka, ipiv,  inv_trans_f, ka, &solve_info);
	pfdaFree(ipiv);
	if(!solve_info){
		pfda_debug_cdl(debugnum_dual_m5_2){pfda_debug_msg("inv_trans_f:\n");printmat(inv_trans_f, *ka, *ka);fflush(stdout);}
		double *tmp = pfdaAlloc_d(*kb**ka, &dp);
		dgemm_(&NoTrans, &NoTrans, kb,ka,ka, &dOne, Lambda , kb, inv_trans_f, ka, &dzero, tmp   ,kb);
		pfdaFree(inv_trans_f);
		pfda_debug_cdl(debugnum_dual_m5_2){pfda_debug_msg("tmp:\n"); printmat(tmp, *kb, *ka);fflush(stdout);}
		dgemm_(&NoTrans, &NoTrans, kb,ka,kb, &dOne, trans_g, kb, tmp        , kb, &dzero, Lambda, kb);
		pfdaFree(tmp);
	} else {
		pfdaFree(inv_trans_f);
		pfda_error("PFDA ERR - pfdaDual_m5_2: could not invert trans_f");
	}
	pfda_debug_cdl(debugnum_dual_m5_2){pfda_debug_msg("leaving pfdaDual_m5_2\n");fflush(stdout);}
}

void pfdaDual_m5(
	int    const * const N,
	int    const * const p,
	double       * const Lambda,
	double       * const alpha,
	double const * const Sigma_aa,
	double       * const Theta_f,
	double       * const Da,
	int    const * const ka,
	double       * const Beta,
	double const * const Sigma_bb,
	double       * const Theta_g,
	double       * const Db,
	int    const * const kb,
	double const * const minimum_variance,
	int const*const dl, double * dp, int*ip)
/* DESCRIPTION
INPUTS:
	dl - integer array for printing debugging information.
WORKING SPACE:
	work_dd - k^2 double
	working space for finding principle components
OUTPUT:
	delta - updated with transformation
	theta - updated with new vectors
	Transformation (k x k) - Q^TTheta
	\[ %tex description
		\Sigma_\delta=\frac{1}{n} \sum\limits_{i=1}{n} \widehat{\delta_i\delta_i\trans}
		\Theta\Sigma_\delta\Theta\trans = Q S_\delta Q\trans
		\delta_i = Q\trans\Theta\delta_i
		Transformation = Q^T\Theta
	\]
MEMORY:
	dp	length = ka^2 + kb^2 + max(
			dpl_5_1  ,
			dpl_5_2 
			)
			
			max(ka,kb) + 2*p^2 + p +  max(  ka, kb, 8 ) * p  (req. for m5_1 which is strictly greater that the req. for m5_2[ ka^2 + kb*ka < 2*p^2])
	ip	length = 6*p
*/
{//code
	pfda_debug_cdl(debugnum_dual_m5){	pfda_debug_msg("pfdaDual_m5 - \n");fflush(stdout);	}
	double * trans_f     = pfdaAlloc_d(*ka**ka,&dp);
	double * trans_g     = pfdaAlloc_d(*kb**kb,&dp);
	

	pfda_debug_cdl(debugnum_dual_m5){	pfda_debug_msg("sub steps\n");fflush(stdout);	}
	pfda_m5_1(N,alpha,ka,Sigma_aa,Theta_f,Da,p,trans_f,minimum_variance,dl,dp,ip);
	pfda_m5_1(N,Beta ,kb,Sigma_bb,Theta_g,Db,p,trans_g,minimum_variance,dl,dp,ip);
	pfdaDual_m5_2(Lambda, trans_f, ka, trans_g, kb, dl,dp,ip);

	pfda_debug_cdl(debugnum_dual_m5){	pfda_debug_msg("freeing memory\n");fflush(stdout);	}
	pfdaFree(trans_f);
	pfdaFree(trans_g);
	pfda_debug_cdl(debugnum_dual_m5){	pfda_debug_msg("leaving pfdaDual_m5\n");fflush(stdout);	}
}

void pfdaDual(
	double * const y, double * const z,
	const double * const B,
	double * const tm,	double * const tn,
	double * const tf,	double * const tg,
	double * const alpha,	double * const beta,
	double * const lambda,
	double * const Da,	double * const Db,
	double * const seps,	double * const sxi,
	double * const Sigma_aa,	double * const Sigma_ab,
	double * const Sigma_bb,
	const int * const nobs,
	const int * const M,
	const int * const N,
	const int * const ka, const int * const kb,
	const int * const p,
	const double * const penalties,
	const double * const ISD,
	const double * const minimum_variance,
	int * const maxI,
	double * const tol,
	const int * const dl,
	double * dp,	int * ip
	)
/*DESCRIPTION:
		This function performs the bivariate fit for paired functional data analysis by Zhou, et al.
	INPUTS:
		y		Mx1 vector of y values   (Also an output for residuals)
		z		Mx1 vector of zvalues  (Also an output for residuals)
		nobs		Nx1 vector of number of observations per subject
		M		number of observations
		N		number of subjects
		ka		number of principle components used for modeling y (ka <= p)
		kb		number of principle components used for modeling z (kb <= p)
		B		Mxp matrix of timepoints evaluated at the orthogonalized basis functions
		p		number of orthogonal basis functions
		minimum_variance		for ridge regression type computations
		--not in find inits--
		penaties	4 vector of penaltie values in the form 	lambda_mu, lambda_nu, lambda_f, lambda_g
		ISD		pxp matrix of the integrated square second derivative of the basis functions
	OUTPUTS:	name|Description(name in R)
		tm		px1 vector estimating mean structore for y (theta_mu)
		tn		px1 vector estimating mean structore for z (theta_nu)
		tf		pxka matrix representing the subject specific portion of the estimate for y(Theta_f)
		tg		pxkb matrix representing the subject specific portion of the estimate for z (Theta_g)
		Da		ka vector representing the variances of alpha (D_alpha)
		Db		kb vector representing the variance of Beta (D_Beta)
		Lambda	kbxka matrix representing the covariance of A and B (Lambda)
		seps		variance of y residuals (sigma_epsilon)
		sxi		variacne of z residuals (sigma_xi)
		alpha		Nxka matrix of PC coefficients of y (alpha)
		Beta		Nxkb matrix of PC coefficients of z(Beta)
		Sigma_aa		kaxkaxN 3 dimentional array of variance components (Sigma_alpha_alpha)
		Sigma_ab		kaxkaxN 3 dimentional array of variance components (Sigma_alpha_beta)
		Sigma_bb		kaxkaxN 3 dimentional array of variance components (Sigma_beta_beta)
	CONTROL VALUES:
		tol		Tolerance for determining convergence
		maxI	on Input: the Maximum number of iterations.
				on output: THe numer of itterations used.
		dl	controls printing of debug information
		incInits	Boolean: does the passing in information include initial values? if true findInits is not called
	MEMORY:
		dp	length = N*p^2+ (ka+kb+2) + ka*kb + ka + kb +
			MAX:
				1. step 1                  	= M + M*max(ka,kb) + 2*max(ka,kb)^2   (exceeded by E)
				2. step 2                 	= p^2+ M+ M*max(ka,kb)
				3. step 3                 	= M + p^2 + p + max(ka)^2 
				4. pfdaDual_m4    	= ka^2 + ka*kb   ( Exceeded by 3)
				5. pfdaDual_m5    	= max(ka,kb) + 2*p^2 + p +  max(  ka, kb, 8 ) * p
				E. pfdaDual_e       	= 2*M+M*ka+M*kb+ka^2+ka*kb+kb^2 + max( 2*max(ni)*max(a,b),4*max(ka,kb)^2)
				(mm4 < m5)
		ip	length = 8 p
*/
{  // Code
	///* declarations //
	//int pcount=0;
	const int kab=*ka**kb, ka2=*ka**ka, kb2=*kb**kb;	//  common constants in computations
	///* Convergence  */
	int I=0; 														//  Itteration number
	double convergenceCriteria=*tol+1;								//  initial value chosen to allow entry into while loop.
	double 	sigma_epsilon_old, sigma_xi_old;						//  Old variables for holding previous state
	double * tm_old     =pfdaAlloc_d(*p, &dp);
	double * tn_old     =pfdaAlloc_d(*p, &dp);
	double * tf_old     =pfdaAlloc_d(*p**ka, &dp);
	double * tg_old     =pfdaAlloc_d(*p**kb, &dp);
	double * Lambda_old =pfdaAlloc_d(kab, &dp);
	double * Da_old     =pfdaAlloc_d(*ka, &dp);
	double * Db_old     =pfdaAlloc_d(*kb, &dp);

	if(*dl){
		pfda_debug_cdl(debugnum_dual_steps)pfda_debug_msg("Entering pfdaDual\n");
		pfda_debug_cdl(debugnum_dual_inputs){
			pfda_debug_msg("INPUTS\n----------------------------------------\n");
			pfda_debug_msg("y:\n"); printmat(y, 1, *M);
			pfda_debug_msg("z:\n"); printmat(z, 1, *M);
			pfda_debug_msg("M:\n%d\n\n",*M);
			pfda_debug_msg("N:\n%d\n\n",*N );
			pfda_debug_msg("ka:\n%d\n\n",*ka);
			pfda_debug_msg("kb:\n%d\n\n",*kb);
			pfda_debug_msg("p:\n%d\n\n",*p);
			pfda_debug_msg("minimum_variance:\n%g\n\n",*minimum_variance);
			pfda_debug_msg("Vector Inputs:\n");
			pfda_debug_msg("nobs:\n"); printmati(nobs, one, *N);//for(i=0;i<*N;i++)pfda_debug_msg("%d\t",nobs[i]); pfda_debug_msg("\n\n");
			pfda_debug_msg("Matrix Arguments:\n");
			pfda_debug_msg("ISD:\n"); printmat(ISD, *p, *p);
			pfda_debug_msg("State Values:\n");
			pfda_debug_msg("tm:\n"); printmat(tm, 1, *p);
			pfda_debug_msg("tn:\n"); printmat(tn, 1, *p);
			pfda_debug_msg("tf:\n"); printmat(tf, *p, *ka);
			pfda_debug_msg("tg:\n"); printmat(tg, *p, *kb);
			pfda_debug_msg("Da:\n"); printmat(Da, 1, *ka);
			pfda_debug_msg("Db:\n"); printmat(Db, 1, *kb);
			pfda_debug_msg("Lambda:\n"); printmat(lambda, *kb, *ka);
			pfda_debug_msg("seps:\n%g\n\n", *seps);
			pfda_debug_msg("sxi:\n%g\n\n", *sxi);
			pfda_debug_msg("alpha:\n"); printmat(alpha, *N, *ka);
			pfda_debug_msg("beta:\n"); printmat(beta, *N, *kb);
			pfda_debug_msg("Criteria:\t%g\nTolerance:\t%g\nMaximun Iterations:\t%d\nI:\t%d\n", convergenceCriteria, *tol, *maxI, I);
		}
		pfda_debug_cdl(debugnum_dual_inputs_large){
			pfda_debug_msg("Sigma_aa:\n"); printmat(Sigma_aa, *N, ka2);
			pfda_debug_msg("Sigma_ab:\n"); printmat(Sigma_ab, *N, kab);
			pfda_debug_msg("Sigma_bb:\n"); printmat(Sigma_bb, *N, kb2);
			pfda_debug_msg("B:\n"); printmat(B, *M, *p);
			}
		fflush(stdout);
	}
	pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Finding Initial values.\n");fflush(stdout);}
	pfdafindinits( y,  z,  nobs, M,  N,  ka, kb, B,  p, minimum_variance,
	   						tm, tn, tf,  tg,  Da,  Db,  lambda,  seps,  sxi,  alpha,  beta,  Sigma_aa,  Sigma_ab,  Sigma_bb); //return values
	setMinVar(seps, &one, minimum_variance);
	setMinVar(sxi, &one, minimum_variance);
	setMinVar(Da, ka, minimum_variance);
	setMinVar(Db, kb, minimum_variance);

	pfda_debug_cdl(debugnum_dual){pfda_debug_msg("Allocating btb(dp=%p)\n",dp);fflush(stdout);}
	double * btb = pfdaAlloc_d(*p**p**N,&dp);
	pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Computing btb\n");fflush(stdout);}
	pfda_computebtb(btb,N,B,M,p,nobs,dl);
	
	pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Entering Convergence Loop\n");fflush(stdout);}
	while(I<*maxI)/*Limited loop with exit condition at end*/{
		pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("entering loop:%d\n", I);fflush(stdout);}

		{ /// * setup convergence */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Setup for convergence\n");fflush(stdout);}
			sigma_epsilon_old = *seps; //teporarily store values for use with convergence
			sigma_xi_old = *sxi;
			dcopy_(&kab, lambda, &one, Lambda_old, &one);
			dcopy_(p, tm, &one, tm_old, &one);
			dcopy_(p, tn, &one, tn_old, &one);
			{int itmp=*p**ka;dcopy_(&itmp, tf, &one, tf_old, &one);}
			{int itmp=*p**kb;dcopy_(&itmp, tg, &one, tg_old, &one);}
			dcopy_(&kab, lambda, &one, Lambda_old, &one);
			dcopy_(ka, Da, &one, Da_old, &one);
			dcopy_(kb, Db, &one, Db_old, &one);
		}
		{ ///* steps 1 */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("steps 1\n");fflush(stdout);}
			pfda_m1(seps, y, nobs,M, N, ka, B, p,minimum_variance, tm, tf, alpha, Sigma_aa, dl, dp	);
			pfda_m1(sxi , z, nobs,M, N, kb, B, p,minimum_variance, tn, tg, beta , Sigma_bb, dl, dp	);
			pfda_debug_cdl(debugnum_dual_estimates){
				pfda_debug_msg("New Parameters\n");
				pfda_debug_msg("seps:\n%# .8g\n\n",*seps);
				pfda_debug_msg("sxi:\n%# .8g\n\n",*sxi);
				fflush(stdout);
			}
		}
		{ ///* steps 2 */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("steps 2\n");fflush(stdout);}
			pfda_m2(tm, y, nobs,M, N, ka, B, p, penalties+0, ISD,tf, seps, alpha, dl, dp);
			pfda_m2(tn, z, nobs,M, N, kb, B, p, penalties+1, ISD,tg, sxi , beta , dl, dp);
			pfda_debug_cdl(debugnum_dual_estimates){
				pfda_debug_msg("New Parameters\n");
				pfda_debug_msg("tm:\n");printmat(tm,*p,one);
				pfda_debug_msg("tn:\n");printmat(tn,*p,one);
				fflush(stdout);
			}
		}
		{ ///* steps 3 */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("steps 3\n");fflush(stdout);}
			pfda_m3(tf, y, nobs, M, N, ka, B, p, penalties+2, ISD,tm, seps, alpha, Sigma_aa, btb, dl, dp, ip);
			pfda_m3(tg, z, nobs, M, N, kb, B, p, penalties+3, ISD,tn, sxi , beta , Sigma_bb, btb, dl, dp, ip);
			pfda_debug_cdl(debugnum_dual_estimates){
				pfda_debug_msg("New Parameters\n");
				pfda_debug_msg("tf:\n");printmat(tf,*p,*ka);
				pfda_debug_msg("tg:\n");printmat(tg,*p,*kb);
				fflush(stdout);
			}
		}
		{ ///* 4 update lambda */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("step 4\n");fflush(stdout);}
			pfdaDual_m4(N,lambda,alpha, ka, beta, kb, Sigma_aa, Sigma_ab,dl, dp);
			pfda_debug_cdl(debugnum_dual_estimates){
				pfda_debug_msg("New Parameters\n");
				pfda_debug_msg("lambda:\n");printmat(lambda,*kb,*ka);
				fflush(stdout);
			}
		}
		{ ///* 5. Orthogonalize  find Da, Db, Update lambda*/
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("step 5\n");fflush(stdout);}
			pfdaDual_m5(N, p, lambda, alpha, Sigma_aa, tf, Da, ka, beta, Sigma_bb, tg, Db, kb, minimum_variance,dl, dp, ip);
			pfda_debug_cdl(debugnum_dual_estimates){
				pfda_debug_msg("New Parameters\n");
				pfda_debug_msg("lambda:\n");printmat(lambda,*kb,*ka);
				pfda_debug_msg("Da:\n");printmat(Da,one,*ka);
				pfda_debug_msg("Db:\n");printmat(Db,one,*kb);
				pfda_debug_msg("tf:\n");printmat(tf,*p,*ka);
				pfda_debug_msg("tg:\n");printmat(tg,*p,*kb);
				fflush(stdout);
			}
		}
		{ ///* E Step */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("step E\n");fflush(stdout);}
			pfdaDual_e(
				y, z, nobs, M, N, ka, kb, B, p,
				tm, tn, tf, tg, Da, Db, lambda,seps, sxi,
				alpha, beta, Sigma_aa, Sigma_ab, Sigma_bb,
				dl, dp ,ip);
			pfda_debug_cdl(debugnum_dual_estimates){
				pfda_debug_msg("New Parameters\n");
				pfda_debug_msg("Sigma_aa:\n");printmat(Sigma_aa,*N,*ka**ka);
				pfda_debug_msg("Sigma_bb:\n");printmat(Sigma_bb,*N,*ka**ka);
				pfda_debug_msg("Sigma_ab:\n");printmat(Sigma_ab,*N,*ka**kb);
				pfda_debug_msg("alpha:\n");   printmat(alpha   ,*N,*ka);
				pfda_debug_msg("beta:\n");    printmat(beta    ,*N,*kb);
				fflush(stdout);
			}
		}
		{ ///* Compute Convergence Criteria */
			pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Computing convergence criteria\n");fflush(stdout);}
			convergenceCriteria = fabs(sigma_epsilon_old-*seps)/MMAX(*seps,*minimum_variance) + fabs(sigma_xi_old-*sxi)/MMAX(*sxi,*minimum_variance);

			pfda_debug_cdl(debugnum_dual_criteria_components){pfda_debug_msg("coputing criteria for seps & sxi:\n");fflush(stdout);}
			daxpy_(p, &mOne, tm, &one, tm_old, &one);
			daxpy_(p, &mOne, tn, &one, tn_old, &one);
			double cctm=0;
			double cctn=0;
			for(int i=0;i<*p;i++){
				cctm+=fabs(tm_old[i])/MMAX(fabs(tm[i]),*minimum_variance);
				cctn+=fabs(tn_old[i])/MMAX(fabs(tn[i]),*minimum_variance);
			}
			convergenceCriteria+=cctm/(*p)+cctn/(*p);

			pfda_debug_cdl(debugnum_dual_criteria_components){pfda_debug_msg("Convergence criteria for tf:\n");fflush(stdout);}
			double cctf=0;
			{int itmp=*p**ka;daxpy_(&itmp, &mOne, tf, &one, tf_old, &one);
			for(int i=0;i<itmp;i++)cctf+=fabs(tf_old[i])/MMAX(fabs(tf[i]),*minimum_variance);
			convergenceCriteria+=cctf/itmp;}

			pfda_debug_cdl(debugnum_dual_criteria_components){pfda_debug_msg("Convergence criteria for tg:\n");fflush(stdout);}
			double cctg=0;
			{int itmp=*p**kb;daxpy_(&itmp, &mOne, tg, &one, tg_old, &one);
			for(int i=0;i<itmp;i++)cctg+=fabs(tg_old[i])/MMAX(fabs(tg[i]),*minimum_variance);
			convergenceCriteria+=cctg/itmp;}

			pfda_debug_cdl(debugnum_dual_criteria_components){pfda_debug_msg("Convergence criteria for lambda:\n");fflush(stdout);}
			daxpy_(&kab, &mOne, lambda, &one, Lambda_old, &one);
			double ccL=0;
			for(int i=0;i<kab;i++)ccL+=fabs(Lambda_old[i])/MMAX(fabs(lambda[i]),*minimum_variance);
			convergenceCriteria+=ccL;

			pfda_debug_cdl(debugnum_dual_criteria_components){pfda_debug_msg("Convergence criteria for Da:\n");fflush(stdout);}
			daxpy_(ka, &mOne, Da, &one, Da_old, &one);
			double ccDa=0;
			for(int i=0;i<*ka;i++)ccDa+=convergenceCriteria+=fabs(Da_old[i])/MMAX(fabs(Da[i]),*minimum_variance);
			convergenceCriteria+=ccDa;

			pfda_debug_cdl(debugnum_dual_criteria_components){pfda_debug_msg("Convergence criteria for Db:\n");fflush(stdout);}
			daxpy_(kb, &mOne, Db, &one, Db_old, &one);
			double ccDb=0;
			for(int i=0;i<*kb;i++)ccDb+=fabs(Db_old[i])/MMAX(fabs(Db[i]),*minimum_variance);
			convergenceCriteria+=ccDb;

			convergenceCriteria/=9;
			///* Debug */
			pfda_debug_cdl(debugnum_dual_criteria_components){
				pfda_debug_msg("Convergence criteria components:\n");
				pfda_debug_msg("seps: \t%5.5g\n", fabs(sigma_epsilon_old-*seps));
				pfda_debug_msg("sxi:  \t%5.5g\n", fabs(sigma_xi_old-*sxi));
				pfda_debug_msg("tm:   \t%5.5g\n", cctm);
				pfda_debug_msg("tn:   \t%5.5g\n", cctn);
				pfda_debug_msg("Tf:   \t%5.5g\n", cctf);
				pfda_debug_msg("Tg:   \t%5.5g\n", cctg);
				pfda_debug_msg("lambda:    \t%5.5g\n", ccL);
				pfda_debug_msg("Da:   \t%5.5g\n", ccDa);
				pfda_debug_msg("Db:   \t%5.5g\n", ccDb);
				fflush(stdout);
			}
			pfda_debug_cdl(debugnum_dual_criteria){pfda_debug_msg("Criteria:%g\n", convergenceCriteria);fflush(stdout);}
		}
		I++;
		if(convergenceCriteria < *tol){
			pfda_debug_cdl(debugnum_dual){pfda_debug_msg("Criteria met leaving loop.%g\n");fflush(stdout);}
			break;
		}
	}
	///* Allow Program to exit normally and leave ensuring convergence to R. */
	{ ///*  Finishing */
		pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Finishing off.\n");fflush(stdout);}
		if(*maxI<=I){pfda_error("EM-algorithm did not converge");}
		pfda_computeResid(y, y,nobs,M, N, ka, B, p,tm, tf, alpha, dl, dp);
		pfda_computeResid(z, z,nobs,M, N, kb, B, p,tn, tg, beta , dl, dp);
		*tol = convergenceCriteria;
		*maxI = I;
	}
	pfda_debug_cdl(debugnum_dual_steps){pfda_debug_msg("Leaving pfdaDual\n");fflush(stdout);}
}

