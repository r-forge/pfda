/*! bin_single.c
Handles the single version of the binary single functional data analysis.
*/
#include "pfda.h"
#include "pfda_gen.h"
#include "pfda_single.h"
#include "pfda_bin_single.h"

// #define 	burninlength 100
// #define		burningenerate 100
// #define		weightedgenerate 10


/*! \defgroup binary_single Univariate Binary

Functions for the univariate binary version of the pfda algorithm
*/

/*! Computes the E step for the single binary version of the pfda algorithm.
\ingroup binary_single
updates as following:
	\f{eqnarray*}{
		\hat\alpha_i &=& \Sigma_{i,\alpha} \Theta_f^T  B_i^T  \left\{E(W_i \vert Y_i)-B_i\theta_\mu\right\} \\
		\Sigma_{i,\alpha} &=& \left( D_\alpha^{-1}+ \theta_f^T B_i^T B_i\theta_f \right)^{-1} \\
		\widehat{\alpha_i\alpha_i^T } 
			&=& \Sigma_{i,\alpha}  + S_{i1}S_{i1}^T + \Sigma_{i,\alpha} \Theta_f^T  B_i^T  E(W_i \vert Y_i) S_{i1}^T  \\
			& & \phantom{\Sigma_{i,\alpha}} + S_{i1} \left\{\Sigma_{i,\alpha} \Theta_f^T  B_i^T  E(W_i \vert Y_i)\right\}^T \\
			& & \phantom{\Sigma_{i,\alpha}} + \Sigma_{i,\alpha} \Theta_f^T  B_i^T
				E(W_i W_i^T \vert Y_i)(\Sigma_{i,\alpha} \Theta_f^T  B_i^T)^T
	\f}
	for each i, with \f$ S_{i1} = -\Sigma_{i,\alpha} \Theta_f^T  B_i^T  B_i\theta_\mu \f$
	
\see pfdaSingle_e()
@MEMORY
	- dp length =	M + M*k + N*k + 2* k^2 + 2 * k + k *max(nobs)
	- ip length  = k
*/
void pfda_bin_single_e(
	double       * const aa_hats,
	double       * const alpha,
	double       * const Sigma_aa,
	double const * const w,
	double const * const ww,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const k,
	double const * const B,
	int    const * const p,
	double const * const tm,
	double const * const tf,
	double const * const Da,
	int * const dl, double *dp , int * ip)
{
	pfda_debug_cdl(debugnum_pfda_bin_single_e){pfda_debug_msg("Entering pfda_bin_single_e\n");fflush(stdout);}
	pfdaSingle_e( alpha, Sigma_aa, w, nobs, M, N, k, B, p, tm, tf, Da, &dOne, dl, dp, ip);
	double * phi = pfdaAlloc_d(*M**k,&dp);
	dgemm_(&NoTrans, &NoTrans, M, k, p, &dOne, B, M, tf, p, &dzero, phi, M);	
	double * omega = pfdaAlloc_d(*M,&dp);
	dgemv_(&NoTrans, M, p, &dOne, B, M, tm, &one, &dzero, omega, &one);
	
	double * S1 = pfdaAlloc_d(*N**k,&dp);
	double * tmp = pfdaAlloc_d(*k,&dp);
	double const * phii=phi, *omegai=omega, * wi = w, *wwi=ww;
	for(int obs=0;obs<*N;obs++){
		pfda_debug_cdl(debugnum_pfda_bin_single_e){pfda_debug_msg("estimating hat{alpha alpha^T} for obs = %d\n",obs);fflush(stdout);}
		dgemv_(&Trans, nobs+obs, k, &dOne, phii, M, omegai, &one, &dzero, tmp, &one);
		dsymv_(&Upper, k, &mOne, Sigma_aa+obs**k**k, k, tmp, &one, &dzero, S1+obs**k, &one);
		pfda_debug_cdl(debugnum_pfda_bin_single_e){pfda_debug_msg("S1[[%d]]\n",obs);printmat(S1+obs**k,*k,one);fflush(stdout);}

		pfda_fillsym(Sigma_aa+obs**k**k,k,dl);
		pfda_bin_cond_aa(aa_hats+obs**k**k, wi, wwi, nobs+obs, M, S1+obs**k, phii, Sigma_aa+obs**k**k, k, dl, dp);
		
		phii += nobs[obs]; omegai += nobs[obs]; wi += nobs[obs]; wwi += nobs[obs]*nobs[obs];
	}	
	pfda_debug_cdl(debugnum_pfda_bin_single_e){pfda_debug_msg("leaving pfda_bin_single_e\n");fflush(stdout);}
}

/*! computes steps 4 and 5 of the binary versions of the pfda algorithm
\MEMORY
	- dp length = k + 2*k^2 + 2*p^2 + p + p*max(k,8)
\ingroup binary_single
\see pfda_m5_0()
*/
void pfda_bin_m5(
	double       * const delta,
	double const * const dd_hats,
	int    const * const N,
	int    const * const k,
	double const * const Sigma_dd,
	double       * const Theta,
	double       * const D,
	int    const * const p,
	double const * const minimum_variance,
	int const * const dl, double * dp, int * ip)
{
	pfda_debug_cdl(debugnum_pfda_bin_m5){ pfda_debug_msg("entering pfda_bin_m5\n"); fflush(stdout);}
	int const k2=*k**k;
	double * sum_dd_hat = pfdaAlloc_d(*k**k, &dp);
	for(int subject=0;subject<*N;subject++){
		pfda_debug_cdl(debugnum_pfda_bin_m5){ pfda_debug_msg("dd_hat[%d]:\n",subject);printmat(dd_hats+subject*k2,*k,*k); fflush(stdout);}
		daxpy_( &k2, &dOne, dd_hats+subject*k2, &one, sum_dd_hat, &one);
	}
	pfda_debug_cdl(debugnum_pfda_bin_m5){ pfda_debug_msg("sum_dd_hat:\n");printmat(sum_dd_hat,*k,*k); fflush(stdout);}
	double * Transformation = pfdaAlloc_d(*k**k,&dp);
	pfda_m5_0( delta, sum_dd_hat, N, k, Sigma_dd, Theta, D, p, Transformation, minimum_variance, dl, dp, ip );
	pfda_debug_cdl(debugnum_pfda_bin_m5){ pfda_debug_msg("leaving pfda_bin_m5\n"); fflush(stdout);}
}


/*! 	@brief This function performs the univariate fit  for functional data analysis using functional principle components
	\ingroup binary_single

	@param	y		[input] Mx1 vector of binary y values
	@param	nobs		[input] Nx1 vector of number of observations per subject
	@param	M		[input] number of observations
	@param	N		[input] number of subjects
	@param	k		[input] number of principle components used for modeling y (k <= p)
	@param	B		[input] Mxp matrix of timepoints evaluated at the orthogonalized basis functions
	@param	p		[input] number of orthogonal basis functions
	@param	minimum_variance		for ridge regression type computations
	@param	lm		[input] lambda_mu, the penalty mutliplier for the mean curve.
	@param	lf		[input] lambda_f, the penalty muultiplier for the principle component curves
	@param	ISD		[input] pxp matrix of the integrated square second derivative of the basis functions
	# OUTPUTS:	name|Description(name in R)
	@param	tm		[output] px1 vector estimating mean structore for y (theta_mu)
	@param	tf		[output] pxka matrix representing the subject specific portion of the estimate for y(Theta_f)
	@param	Da		[output] k vector representing the variances of alpha (D_alpha)
	@param	alpha		[output] Nxka matrix of PC coefficients of y (alpha)
	@param	Sigma_aa	[output] kaxkaxN 3 dimentional array of variance components (Sigma_alpha_alpha)
	CONTROL VALUES:
	@param	convergencetol		Tolerance for determining convergence of overall algorithm
	@param	maxI	on Input: the Maximum number of iterations.
				          on output: THe numer of itterations used.
	@param 	burninlength	The length of the burnin period 
	@param burningenerate	The number of replicates to generate during the burnin period
	@param weightedgenerate  The number of replicates to generate after the burnin period

	@param	dl	controls printing of debug information
	@param	dp	pool of preallocated double memory.
	@param	ip	pool of preallocated integer memory.
	@MEMORY
		- dp	length = M + sum(nobs^2) + N*k^2 + N*p^2 + p + p*k + k
			\n + MAX:
				- inits                      = 2*p^2 + M + p*N + 8*p
				- step W                	= ni*kr + kr + (ni * k + ni + k^2 + k + p + p^2 + k*p + 3 k)
				- step 1/E              	= M + M*k + N*k + 2* k^2 + 2 * k + k *ni
				- step 2                 	= p^2+ M+ M*k
				- step 3                 	= M + p^2 + p + k^2 
				- step 4                  	= k + 2*k^2 + 2*p^2 + p + p*max(k,8)
				.
				where ni = max(nobs)
		- ip	length = 8 p
	\callgraph 
*/
void pfda_bin_single(
	/* input Values */
	int    const * const y,
	int    const * const  nobs,
	int    const * const  M,
	int    const * const  N,
	int    const * const  k,
	double const * const B,
	int    const * const  p,
	double const * const lm,
	double const * const lf,
	double const * const ISD,
	/* State Values */
	double       * const tm,
	double       * const tf,
	double       * const Da,
	double       * const alpha,
	double       * const Sigma_aa,
	/* Control Values */
	double const * const minimum_variance, 
	double       * const convergencetol,
	int * const maxI,
	int const * const burninlength,
	int const * const burningenerate,
	int const * const weightedgenerate,
	int * const dl,
	double *dp , int * ip)
{ pfda_debug_step;
if(*dl){
	pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Entering pfda_bin_single\n");fflush(stdout);}
	pfda_debug_cdl(debugnum_bin_single_inputs){
		pfda_debug_msg("INPUTS\n----------------------------------------\n");
		pfda_debug_msg("M:\t%d\n",*M);
		pfda_debug_msg("N:\t%d\n",*N );
		pfda_debug_msg("k:\t%d\n",*k);
		pfda_debug_msg("p:\t%d\n",*p);
		pfda_debug_msg("lm:\t%g\n",*lm);
		pfda_debug_msg("lf:\t%g\n",*lf);
		pfda_debug_msg("minimum_variance:\n%g\n",*minimum_variance);
		pfda_debug_msg("Convergence Tolerance:\t%g\n",  *convergencetol);
		pfda_debug_msg("Maximun Iterations:\t%d\n",  *maxI);
		pfda_debug_msg("burninlength:\t%d\n",*burninlength);
		pfda_debug_msg("burningenerate:\t%d\n",*burningenerate);
		pfda_debug_msg("weightedgenerate:\t%d\n",*weightedgenerate);
		fflush(stdout);

		pfda_debug_msg("Vector Arguments:\n");
		pfda_debug_msg("y:\n"); printyveci(y, nobs, N);
		pfda_debug_msg("nobs:\n"); for(int i=0;i<*N;i++)pfda_debug_msg("%d\t",nobs[i]); pfda_debug_msg("\n\n");
		fflush(stdout);

		pfda_debug_msg("Outputs:");
		pfda_debug_msg("tm:\n"); printmat(tm, 1, *p);
		pfda_debug_msg("tf:\n"); printmat(tf, *p, *k);
		pfda_debug_msg("Da:\n"); printmat(Da, 1, *k);
		pfda_debug_msg("alpha:\n"); printmat(alpha, *N, *k);
		fflush(stdout);
	}
	pfda_debug_cdl(debugnum_bin_single_inputs_large){
		pfda_debug_msg("Large Inputs:");
		pfda_debug_msg("ISD:\n"); printmat(ISD, *p, *p);
		pfda_debug_msg("B:\n"); printmat(B, *M, *p);
		pfda_debug_msg("Sigma_aa:\n"); printmat(Sigma_aa, *N, *k**k);
	}
}
GetRNGstate();
double * w = pfdaAlloc_d(*M,&dp);
int size_ww = 0;for(int i=0;i<*N;i++)size_ww+=nobs[i]*nobs[i];
double * ww = pfdaAlloc_d(size_ww, &dp);
double * aa_hats = pfdaAlloc_d(*N**k**k,&dp);//used to hold the conditional expectations used in steps 3 and step 5
double * btb = pfdaAlloc_d(*p**p**N,&dp);
pfda_computebtb(btb,N,B,M,p,nobs,dl);

{ // find initial values
	pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Estimating initial values\n");fflush(stdout);}
	for(int i=0;i<*M;i++)w[i]=(double)y[i];
	pfda_s_i(tm, tf, alpha, Da, aa_hats, NULL, w,  nobs, M, N, k,  B,  btb, p, minimum_variance, dl, dp, ip);
	setMinVar(Da, k, minimum_variance);
}

int I=0;
double const * const seps=&dOne;
double   convergenceCriteria=0;
double * tmOld  = pfdaAlloc_d(*p    , &dp);
double * tf_old = pfdaAlloc_d(*p**k, &dp);
double * Da_old = pfdaAlloc_d(*k   , &dp);
while(I < *maxI)/*  */{
	{/* setup for convergence */
		pfda_debug_cdl(debugnum_bin_single_steps){
			pfda_debug_msg("Entering loop: %d\n",I);
			pfda_debug_msg("Setup for convergence\n");
			fflush(stdout);
		}
		dcopy_(p, tm, &one, tmOld, &one);
		int pka=*p**k;dcopy_(&pka, tf, &one, tf_old, &one);
		dcopy_(k, Da, &one, Da_old, &one);
	}
	{/* W Step - W|y / WW|Y */
		pfda_debug_cdl(debugnum_pfda_bin_single_break_before_ew){
			pfda_debug_msg("debuging information has indicated to exit prematurly before the W step.\n");
			pfda_debug_msg("Final States of w and ww\n");
			pfda_debug_msg("w:\n");printyvec(w,nobs,N);
			pfda_debug_msg("ww:\n");double * wwi=ww;for(int i =0; i<*N;i++){
				pfda_fillsym(wwi, nobs+i, dl);
				pfda_debug_msg("[[%d]]\n",i);printmat(wwi,nobs[i],nobs[i]);wwi+=nobs[i]*nobs[i];
			}
			fflush(stdout);
			return;
		}
		pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("W Step\n");fflush(stdout);}
		int const *  kr = (I<*burninlength)?burningenerate:weightedgenerate;
		double weight = (I<*burninlength)?1:(10.0/(10.0+(double)I));
		pfda_bin_single_approximate_moments( w, ww, y, nobs, M, N, B, p, k, tm, tf, Da, &weight, kr, dl, dp, ip);
	}
	{/* step 1/E - alpha / Sigma_aa / aa_hats*/
		pfda_debug_cdl(debugnum_pfda_bin_single_break_before_e){pfda_debug_msg("debuging information has indicated to exit prematurly before the e step.\n");fflush(stdout);return;}
		pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Step E\n");fflush(stdout);}
		pfda_bin_single_e( aa_hats, alpha, Sigma_aa, w, ww, nobs, M, N, k, B, p, tm, tf, Da, dl, dp, ip);
		pfda_debug_cdl(debugnum_bin_single_estimates){
			pfda_debug_msg("alpha: \n");printmat(alpha,*p,*k);
			pfda_debug_msg("Sigma_aa: \n");printmat(Sigma_aa,*N,*k**k);
			for(int i=0;i<*N;i++){pfda_debug_msg("aa_hat[[%d]]: \n",i);printmat(aa_hats+i**k**k,*k,*k);}
			fflush(stdout);}
	}
	{/* step 2 - tm*/
		pfda_debug_cdl(debugnum_pfda_bin_single_break_before_2){pfda_debug_msg("debuging information has indicated to exit prematurly before step 2.\n");fflush(stdout);return;}
		pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Step 2\n");fflush(stdout);}
		pfda_m2(tm,w,nobs,M,N,k,B,p,lm,ISD,tf,seps,alpha,dl, dp);
		pfda_debug_cdl(debugnum_bin_single_estimates){pfda_debug_msg("tm: \n");printmat(tm,*p,one);fflush(stdout);}
	}
	{/* step 3 -tf */
		pfda_debug_cdl(debugnum_pfda_bin_single_break_before_3){pfda_debug_msg("debuging information has indicated to exit prematurly before step 3.\n");fflush(stdout);return;}
		pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Step 3\n");fflush(stdout);}
		pfda_m3_core(tf, aa_hats, w, nobs, M, N, k, B, p, lf, ISD, tm, seps, alpha, btb, dl, dp, ip);
		pfda_debug_cdl(debugnum_bin_single_estimates){pfda_debug_msg("tf: \n");printmat(tf,*p,*k);fflush(stdout);}
	}
	{/* step 4 - Da/ orthogonalize tf*/
		pfda_debug_cdl(debugnum_pfda_bin_single_break_before_4){pfda_debug_msg("debuging information has indicated to exit prematurly before step 4.\n");fflush(stdout);return;}
		pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Step 4&5\n");fflush(stdout);}
		pfda_bin_m5( alpha, aa_hats, N, k, Sigma_aa, tf, Da, p, minimum_variance, dl, dp, ip);
		pfda_debug_cdl(debugnum_bin_single_estimates){
			pfda_debug_msg("tf: \n");printmat(tf,*p,*k);
			pfda_debug_msg("Da: \n");printmat(Da,*k,one);
			pfda_debug_msg("alpha: \n");printmat(alpha,*p,*k);
			fflush(stdout);}
	}
	{/* Compute Convergence Criteria */
		pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Computing convergence Criteria\n");fflush(stdout);}

		daxpy_(p, &mOne, tm, &one, tmOld, &one);
		double cctm=0;
		for(int i=0;i<*p;i++){
			cctm+=fabs(tmOld[i]/tm[i]);
		}
		convergenceCriteria=cctm;

		double cctf=0;
		{int pka=*p**k;daxpy_(&pka, &mOne, tf, &one, tf_old, &one);
		for(int i=0;i<pka;i++)cctf+=fabs(tf_old[i]/tf[i]);
		convergenceCriteria+=cctf;}

		//daxpy_(k, &mOne, Da, &one, Da_old, &one);
		double ccDa=0;
		for(int i=0;i<*k;i++)ccDa+=fabs(Da_old[i]-Da[i])/Da[i];
		convergenceCriteria+=ccDa;

		/* Debug */
		pfda_debug_cdl(debugnum_bin_single_convergence_components){
				pfda_debug_msg("Convergence criteria components:\n");
				pfda_debug_msg("tm:   \t%5.5g\n", cctm);
				pfda_debug_msg("Tf:   \t%5.5g\n", cctf);
				pfda_debug_msg("Da:   \t%5.5g\n", ccDa);
				fflush(stdout);
		}
		pfda_debug_cdl(debugnum_bin_single_convergence){pfda_debug_msg("Criteria:%g\n", convergenceCriteria);fflush(stdout);}
	}
	I++;
	if(convergenceCriteria < *convergencetol)break;
}//end while loop
{/* Finishing */
	pfda_debug_cdl(debugnum_bin_single_steps){pfda_debug_msg("Finishing\n");fflush(stdout);}
	if(*maxI<=I){pfda_error("EM-algorithm did not converge");}
	//pfda_computeResid(y, y, nobs, M, N, k, B, p, tm, tf, alpha, dl, dp);
	*convergencetol = convergenceCriteria;
	*maxI = I;
}
pfda_debug_cdl(debugnum_bin_single_steps){ pfda_debug_msg("Leaving pfdaSingle\n");fflush(stdout); }
PutRNGstate();
}//end pfdaSingle
