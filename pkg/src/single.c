/*! \file single.c

// by Andrew Redd
// Part of paired functional data analysis through principle component analysis.
// This file performs the em algorithm for a single Curve.  Very similar to the paired algorithm and possible should be combined.
*/
#include "pfda_error.h"
#include "pfda_gen.h"
#include "pfda_memory.h"
#include "pfda_single.h"
#include "pfda_utils.h"
#include <R_ext/Lapack.h>

/*! compute the residual vector
	\f[
	y_i = y_i - B_i tm - B_i tf alpha_i
	\f]
@MEMORY
	dp = M*k if tf is not NULL
*/
void pfda_computeResid(
	double       * const Ry, 		//on exit y will hold residuals.
	double const * const y, 		//on exit y will hold residuals.
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const k,
	double const * const B,
	int    const * const p,
	double const * const tm,
	double const * const tf,
	double const * const alpha,
	int const * const dl, double*dp)
{pfda_debug_step
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(tm)dgemv_(&NoTrans, M, p, &mOne, B, M, tm, &one, &dOne, Ry, &one); //y = (y-B*tm)
	if(tf && alpha && nobs && N){
		double * psi = pfdaAlloc_d(*M**k,&dp);
		dgemm_(&NoTrans, &NoTrans, M, k, p, &dOne, B, M, tf, p, &dzero, psi, M);//  psi = B*tg
		double * Ryi=Ry; double * psii=psi;
		for(int obsnum=0; obsnum<*N; obsnum++){
			dgemm_(&NoTrans, &Trans, &(nobs[obsnum]), &one, k, &mOne, psii, M, alpha+obsnum, N, &dOne, Ryi, M);//Ry-=B*Tf*alpha[i]]
			Ryi+=nobs[obsnum]; psii+=nobs[obsnum];
		}
	}
}

/*! Compute residual vector includes possibility for extra Z variables
 *
 */
void single_c_resid(
	double       * const Ry,
	double const * const y,
	double const * const Z,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const kz,
	int    const * const k,
	double const * const B,
	int    const * const p,
	double const * const tz,
	double const * const tm,
	double const * const tf,
	double const * const alpha,
	int const * const dl, double*dp)
{pfda_debug_step
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(Z && kz && tz && *kz)dgemv_(&NoTrans, M, kz, &mOne, Z, M, tz, &one, &dOne, Ry, &one); //y = (y-Z*tz)
	pfda_computeResid( Ry, Ry, nobs, M, N, k, B, p, tm, tf, alpha, dl, dp);
}

/*! takes the matrix B partintions it by nobs and computes the inner product of B_i^t B_i and stores it in btb[i] as upper symetric matrices
*/
void pfda_computebtb(
	double *btb,
	int const * const N,
	double const * const B,
	int const * const M,
	int const * const p,
	int const * const nobs,
	int const * const dl)
{pfda_debug_step
	int subjectnum=0;double const * Bi=B;
	for(subjectnum=0;subjectnum<*N;subjectnum++){
		dsyrk_(&Upper, &Trans, p, nobs+subjectnum, &dOne, Bi, M, &dzero, btb+subjectnum**p**p, p);//btb[subjectnum]=B_i'B_i
		Bi+=nobs[subjectnum];
	}
}

/*!  Step 1 - find sigma_epsilon
@MEMORY
	-dp = M + M*k + 2*k^2
*/
void pfda_m1(
	double       * const sigma,
	/* Inputs Values */
	double const * const y,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const ka,
	double const * const B,
	int    const * const p,
	double const * const minV,
	double const * const tm,
	double const * const tf,
	double const * const alpha,
	double const * const Saa,
	/* Auxilliary */
	int const * const dl, double * dp
	)
{pfda_debug_step
	// double *BTfi2=NULL, *BTfi3=NULL;
	pfda_debug_cdl(debugnum_m1){
		pfda_debug_msg("pfda_m1 - \n");
		pfda_debug_msg("y:\n");printyvec(y,nobs,N);
		pfda_debug_argmati(nobs,1,*N);
		pfda_debug_arg(*M);
		pfda_debug_arg(*N);
		pfda_debug_arg(*ka);
		pfda_debug_arg(*p);
		pfda_debug_arg(*minV);
		pfda_debug_argmat(tm,*p,1);
		pfda_debug_argmat(tf,*p,*ka);
		pfda_debug_argmat(alpha,*N,*ka);
		fflush(stdout);
	}
	pfda_debug_cdl(debugnum_m1_inputs){
		pfda_debug_argmat(B,*M,*p);
		int sn;for(sn=0;sn<*N;sn++){
			pfda_debug_arg(sn);pfda_debug_argmat(Saa+sn**ka**ka,*ka,*ka);
		}
		fflush(stdout);
	}	///* memory allocation */
	int maxn = nobs[1], i;for(i=1;i<*N;i++)maxn=MMAX(maxn,nobs[i]);
	double       * Ry     = pfdaAlloc_d(*M     , &dp);
	double       * phi    = pfdaAlloc_d(*M**ka , &dp);
	double       * phiphi    = pfdaAlloc_d(*ka**ka, &dp);
	double       * Sigmaphiphi    = pfdaAlloc_d(*ka**ka, &dp);
	double       * phii   = phi;
	double       * Ryi    = Ry;
	double const * Sigma_i= Saa;
	int subjectnum;

	pfda_debug_cdl(debugnum_m1){pfda_debug_msg("Compute Ry - \n");fflush(stdout);}
	dcopy_(M,y,&one,Ry,&one);
	dgemv_(&NoTrans, M, p, &mOne, B, M, tm, &one, &dOne, Ry, &one); //Ry -= B*tm  (y-B*tm)
	pfda_debug_cdl(debugnum_m1){pfda_debug_msg("Ry = y-B*tm:\n");printyvec(Ry,nobs,N);}
	dgemm_(&NoTrans, &NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);//compute phi = B*tf
	pfda_debug_cdl(debugnum_m1_inputs){pfda_debug_argmat(phi,*M,*ka);}

	pfda_debug_cdl(debugnum_m1){pfda_debug_msg("Enter subject loop - \n");fflush(stdout);}
	*sigma = 0.0;
	Ryi=Ry; Sigma_i=Saa; phii=phi;
	for(subjectnum=0; subjectnum<*N; subjectnum++){
		pfda_debug_cdl(debugnum_m1_inloop){pfda_debug_arg(subjectnum);fflush(stdout);}
		dgemm_(&NoTrans, &Trans, nobs+subjectnum, &one, ka, &mOne, phii, M, alpha+subjectnum, N, &dOne, Ryi, M);//Ry-=B*Tf*alpha[i]]
		dsyrk_(&Upper, &Trans, ka, nobs+subjectnum, &dOne, phii, M, &dzero, phiphi, ka);
		pfda_debug_cdl(debugnum_m1_inloop){pfda_debug_argmat(phiphi,*ka,*ka);fflush(stdout);}
		dsymm_(&Left, &Upper, ka, ka, &dOne, phiphi, ka, Sigma_i, ka, &dzero, Sigmaphiphi, ka);
		pfda_debug_cdl(debugnum_m1_inloop){pfda_debug_argmat(Sigmaphiphi,*ka,*ka);fflush(stdout);}
		*sigma += tr(Sigmaphiphi, *ka); //add traces
		pfda_debug_cdl(debugnum_m1_inloop){pfda_debug_msg("sigma(partial sum): %g\n", *sigma);}
		// /* Increment partial pointers */
		phii+=nobs[subjectnum];
		Ryi +=nobs[subjectnum];
		Sigma_i+=*ka**ka;
	}
	pfda_debug_cdl(debugnum_m1){
		pfda_debug_msg("Ry=(y-B*tm-B*Tf*alpha):\n");pfda_debug_argmat(Ry,1,*M);
		pfda_debug_msg("sigma(sum of traces): %g\n", *sigma);
	}

	*sigma += ddot_(M, Ry, &one, Ry, &one);//sigma += Ry^t*Ry
	pfda_debug_cdl(debugnum_m1){pfda_debug_msg("sigma(with Ry^t*Ry): %g\n", *sigma);}

	*sigma /= *M;//divide by number of observations
	setMinVar(sigma, &one, minV);
	pfda_debug_cdl(debugnum_m1){
		pfda_debug_msg("Step 1 results:\n");
		pfda_debug_msg("sigma:\t%g\n", *sigma);
	}
}

/*! Step 2 - stimates the mean curve theta_mu
 *
 * step 2 of the M step of the EM algorithm for paired functional data analysis by Zhou, et al.
 * @MEMORY
 * 	-dp = p^2+ M+ M*ka
 */
void pfda_m2(
	double * const tm,
	/* input Values */
	double const * const y,
	int const * const  nobs,
	int const * const  M,
	int const * const  N,
	int const * const  ka,
	double const * const B,
	int const * const  p,
	double const * const lm,
	double const * const K,
	double const * const tf,
	double const * const sigma,
	double const * const alpha,
	int const * const dl, double * dp)
{pfda_debug_step
	pfda_debug_cdl(debugnum_m2){
		pfda_debug_msg("pfda_m2 - \n");
		fflush(stdout);
		pfda_debug_msg("y:\n");printyvec(y,nobs,N);
		pfda_debug_argmati(nobs,one,*N);
		pfda_debug_arg(*M);
		pfda_debug_arg(*N);
		pfda_debug_arg(*ka);
		pfda_debug_arg(*p);
		pfda_debug_arg(*lm);
		pfda_debug_argmat(tf,*p,*ka);
		pfda_debug_arg(*lm);
		pfda_debug_argmat(alpha,*N,*ka);
		fflush(stdout);
	}
	pfda_debug_cdl(debugnum_m2_inputs_large){
		pfda_debug_argmat(B,*M,*p);
		pfda_debug_argmat(K,*p,*p);
		fflush(stdout);
	}
	double * Ry        = pfdaAlloc_d(*M, &dp);
	dcopy_(M, y, &one, Ry, &one); //set Ry=y

	double       * phi    = pfdaAlloc_d(*M**ka ,&dp);
	dgemm_(&NoTrans, &NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);//compute phi = B*tf

	for(int i=0;i<*p;i++){tm[i]=0;}


	double const * Bi=B; double * Ryi=Ry; double * phii=phi;
	int subjectnum=0;
	for(subjectnum=0; subjectnum<*N; subjectnum++){
		dgemv_(&NoTrans, &(nobs[subjectnum]), ka,  &mOne, phii, M, &(alpha[subjectnum]), N, &dOne, Ryi, &one);
		dgemv_(&Trans, &(nobs[subjectnum]), p, &dOne, Bi, M, Ryi, &one, &dOne, tm, &one); //  += B_i'Ryi
		// /* Increment partial pointers */
		Bi    += nobs[subjectnum];
		phii  += nobs[subjectnum];
		Ryi   += nobs[subjectnum];
	}
	pfdaFree(Ry);
	pfdaFree(phi);
	pfda_debug_cdl(debugnum_m2){
		pfda_debug_msg("Ry:\n");printyvec(Ry,nobs,N);
		pfda_debug_argmat(tm, *p, one);
		fflush(stdout);
	}


	pfda_debug_cdl(debugnum_m2){pfda_debug_msg("PING!\n");fflush(stdout);}
	double * BtB  = pfdaAlloc_d(*p**p, &dp);
	pfda_debug_cdl(debugnum_m2){pfda_debug_msg("PING!\n");fflush(stdout);}
	dsyrk_(&Upper, &Trans, p, M, &dOne, B, M, &dzero, BtB, p);    // BtB = B^t B
	pfda_debug_cdl(debugnum_m2){pfda_debug_msg("PING!\n");fflush(stdout);}
	for(int i=0;i<*p;i++)for(int j=0;j<=i;j++){  //BtB += sigma * lm * K
		BtB[i**p+j]+=*sigma**lm*K[i**p+j];
	}
	pfda_debug_cdl(debugnum_m2){pfda_debug_msg("sum B^t B + sigma*lm*K (left side of solve) = \n");pfda_debug_argmat(BtB,*p,*p);fflush(stdout);}

	int sr = 0;
	dposv_(&Upper, p, &one, BtB, p, tm, p, &sr);
	pfdaFree(BtB);
	if(sr)pfda_error("PFDA ERR: Leading minor of order %i is not positive definite (theta_mu)",sr);

	pfda_debug_cdl(debugnum_m2){
		pfda_debug_argmat(tm, *p, one);
		pfda_debug_msg("leaving pfda_m2\n");
		fflush(stdout);
	}
}

/*! most interior of the functions for step 3 of the M step

Contains the most interoir of the 3 functions that comprise the M step of the Algorithm.
adds onto the running totals for left and right.

@MEMORY
	- dp length = p
*/
void pfda_m3_for_subjectnum(
	double       * right,
	double       * left,
	double const * const Ryi,
	int    const * const n_i,
	int    const * const M,
	int    const * const k,
	double const * const Bi,
	double const * const tf,
	int    const * const p,
	double const * const alpha_ij,
	double const * const aa_hat_i,
	double const * const btb_i,
	int    const * const subject_num,
	int    const * const estcol_num,
	int const * const dl,
	double * dp
	)
{	pfda_debug_step
	int const estcolnum = *estcol_num;
	pfda_debug_cdl(debugnum_m3_level2){pfda_debug_arg(*subject_num);fflush(stdout);}
	double * sum_cols = pfdaAlloc_d(*p,&dp);
	pfda_debug_cdl(debugnum_m3_level2){pfda_debug_argmat(aa_hat_i, *k,*k);fflush(stdout);}
	pfda_uppersymadd(p, left, p, aa_hat_i+estcolnum*(*k+1), btb_i, p, dl);
	pfda_debug_cdl(debugnum_m3_level2){pfda_debug_argmat(left, *p,*p);fflush(stdout);}
	for(int cellnum=0;cellnum<*p;cellnum++){sum_cols[cellnum]=0;}
	pfda_debug_cdl(debugnum_m3_level2){pfda_debug_msg("entering non-estimating column loop loop\n");fflush(stdout);}
	for(int colnum=0;colnum<*k;colnum++)if(colnum!=estcolnum){
		daxpy_(p, aa_hat_i+MMAX(estcolnum,colnum)**k+MMIN(estcolnum,colnum), tf+colnum**p, &one, sum_cols, &one);
		pfda_debug_cdl(debugnum_m3_level3){
			pfda_debug_msg("adding column # %d\n",colnum);
			pfda_debug_msg("aa_hat_i offset: %d\n",MMAX(estcolnum,colnum)**k+MMIN(estcolnum,colnum));
			pfda_debug_msg("multiplier: %9.8g\n",*(aa_hat_i+MMAX(estcolnum,colnum)**k+MMIN(estcolnum,colnum)));
			pfda_debug_argmat(tf+colnum**p,*p,one);
			pfda_debug_argmat(sum_cols,*p,one);
			fflush(stdout);}
	}
	pfda_debug_cdl(debugnum_m3_level2){
		pfda_debug_argmat(btb_i,*p,*p);
		pfda_debug_argmat(sum_cols,*p,one);
		fflush(stdout);}
	dsymv_(&Upper, p, &mOne, btb_i, p, sum_cols, &one, &dOne, right, &one);
	pfda_debug_cdl(debugnum_m3_level2){pfda_debug_msg("right(partial btb x sum_cols):\n");pfda_debug_argmat(right,*p,one);fflush(stdout);}
	dgemv_(&Trans, n_i, p, alpha_ij, Bi, M, Ryi, &one, &dOne, right, &one);
	pfda_debug_cdl(debugnum_m3_level2){pfda_debug_argmat(right,*p,one);fflush(stdout);}
}

/*! middle fucntion for the 3 comprising the 3rd step of the M step in the algorithm.

estimates the new column provided by right

@MEMORY
	- dp length = p^2 + p
	- ip length = p
*/
void pfda_m3_for_estcol(
	double * const right,
	/* input Values */
	double const * const aa_hats,
	double const * const Ry,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const k,
	double const * const B,
	int    const * const p,
	double const * const lf,
	double const * const K,
	double const * const tm,
	double const * const tf,
	double const * const sigma,
	double const * const alpha,
	double const * const btb,
	int const * const estcol_num,
	int const * const dl, double * dp, int * ip)
{pfda_debug_step
	int const estcolnum=*estcol_num;
	pfda_debug_cdl(debugnum_m3_level1){pfda_debug_msg("estimating column =%d\n",estcolnum);fflush(stdout);}
	//double * right=tf+estcolnum**p;
	int const p2 = *p**p;
	double * left     = pfdaAlloc_d(p2,&dp);
	for(int cellnum=0;cellnum<*p;cellnum++){right[cellnum]=0.0;}
	for(int cellnum=0;cellnum<p2;cellnum++){left[cellnum]=0.0;}
	pfda_debug_cdl(debugnum_m3_level1){pfda_debug_msg("entering subject loop\n");fflush(stdout);}
	double const *Bi=B, *Ryi=Ry;
	for(int subjectnum=0;subjectnum<*N;subjectnum++){
		pfda_m3_for_subjectnum(right, left, Ryi, nobs+subjectnum, M, k, Bi, tf, p,
			alpha+estcolnum**N+subjectnum, aa_hats+subjectnum**k**k,
			btb+subjectnum*p2, &subjectnum, &estcolnum, dl, dp);
		Bi+=nobs[subjectnum];
		Ryi+=nobs[subjectnum];
	}
	double penalty_multiplier = *sigma**lf;
	pfda_uppersymadd(p, left, p, &penalty_multiplier, K, p, dl);

	pfda_debug_cdl(debugnum_m3_level1){pfda_debug_argmat(left, *p,*p);fflush(stdout);}
	pfda_debug_cdl(debugnum_m3_level1){pfda_debug_argmat(right,*p,one);fflush(stdout);}
	int sr=0;
	int lwork = -1;
	double thiswork =0;
	double * work = &thiswork;
	int * ipiv = pfdaAlloc_i(*p, &ip);
	dsysv_(&Upper, p, &one, left, p, ipiv,  right, p, work, &lwork, &sr);  // work size querry
	lwork = *work;
	work = Calloc(lwork, double);  // since size is unpredictable dynamic allocation must be used.
	dsysv_(&Upper, p, &one, left, p, ipiv,  right, p, work, &lwork, &sr);  // actually sovle the system
	Free(work);
	if(sr<0){
		pfda_error("PFDA ERR (pfda_m3): dsysv return error code: %d on estimating column: %d",sr,estcolnum+1);
	}
}

/*! Step 3 -  Estimate Theta_f
@MEMORY:
	- dp  length =  M + k^2 * N + p^2 + p
	- ip  length = p
 */
void pfda_m3_core(
	double * const tf,
	/* input Values */
	double const * const aa_hats,
	double const * const y,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const k,
	double const * const B,
	int    const * const p,
	double const * const lf,
	double const * const K,
	double const * const tm,
	double const * const sigma,
	double const * const alpha,
	double const * const btb,
	int const * const dl, double * dp, int * ip)
{pfda_debug_step
	pfda_debug_cdl(debugnum_m3){
		pfda_debug_msg("entering pfda_m3_core\n");
		pfda_debug_arg(*lf);
		pfda_debug_argmat(tm,*p,one);
		pfda_debug_argmat(tf,*p,*k);
		pfda_debug_argmat(alpha,*N,*k);
		pfda_debug_argmat(B,*M,*p);
		fflush(stdout);
	}

	double * Ry = pfdaAlloc_d(*M,&dp);

	dcopy_(M, y, &one, Ry, &one); //set Ry=y
	dgemv_(&NoTrans, M, p, &mOne, B, M, tm, &one, &dOne, Ry, &one); //Ry -= B*tm  (y-B*tm)

	pfda_debug_cdl(debugnum_m3){pfda_debug_msg("loop prep.\n");fflush(stdout);}

	pfda_debug_cdl(debugnum_m3){pfda_debug_msg("entering estimating column loop\n");fflush(stdout);}
	for(int estcolnum=0;estcolnum<*k;estcolnum++){
		pfda_m3_for_estcol(tf+estcolnum**p, aa_hats, Ry, nobs, M, N, k, B, p, lf,
			K, tm, tf, sigma, alpha,  btb, &estcolnum, dl, dp, ip);
	}
	pfda_debug_cdl(debugnum_m3){pfda_debug_argmat(tf,*p,*k);fflush(stdout);}
	pfda_debug_cdl(debugnum_m3){pfda_debug_msg("leaving pfda_m3_core\n");fflush(stdout);}
}

/*! @brief wrapper for the pfda_m3_core that does all the work for the single linear alorithm

only adds computing the aa_hats then passes them on.
@MEMORY
	- dp  length =  M + k^2 * N + p^2 + p
	- ip  length = p
*/
void pfda_m3(
	double * const tf,
	/* input Values */
	double const * const y,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const k,
	double const * const B,
	int    const * const p,
	double const * const lf,
	double const * const K,
	double const * const tm,
	double const * const sigma,
	double const * const alpha,
	double const * const Saa,
	double const * const btb,
	int const * const dl, double * dp, int * ip)
{pfda_debug_step
	int const k2 = *k**k;
	double * aa_hats = pfdaAlloc_d(k2**N, &dp);
	for(int subject=0;subject<*N;subject++){
		pfda_cond_dd( aa_hats+subject**k**k, alpha+subject, N, Saa+subject*k2, k, dl);
		pfda_debug_cdl(debugnum_m3){pfda_debug_argmat(aa_hats+subject**k**k,*k,*k);fflush(stdout);}
	}
	pfda_m3_core( tf, aa_hats, y, nobs, M, N, k, B, p, lf, K, tm, sigma, alpha, btb, dl, dp, ip);
}

/*! computes the core of the computation for the m5 step for the single

Among the single version the binary and linear they are differentiatied in separate functions by specifying sum_dd_hat  \f$ =\sum\limits_{i=1}^N \widehat{\delta_i\delta_i^T} \f$
where \f$ \delta \f$ is a place holder for \f$ \alpha \f$ or \f$ \beta \f$.

	@param delta [output]
	@param sum_dd_hat already summed up
	@param N Number of subjects
	@param k
	@param Sigma_dd
	@param Theta [output]
	@param D [output]
	@param p
	@param Transformation  [output] k*k matrix to hold the transformation matrix
	@param minV
	@param dl
	@param dp
	@param ip .

Computes and performs the Eigen value decomposition to orthogonalize the system.
	  \f{eqnarray*}{
		\Sigma_\delta                            &=&\frac{1}{n} \sum\limits_{i=1}{n} \widehat{ \delta_i \delta_i^T  }	\\
		\Theta\Sigma_\delta\Theta^T &=& Q D Q^T							\\
		\delta_i                                       &=& Q^T\Theta\delta_i						\\
		Transformation &=& Q^T\Theta
	\f}
@MEMORY
	- dp length = k + 2*p^2 + 9*p +  k*p ,
	- ip length = 6*p
*/
void pfda_m5_0(
	double       * const delta,
	double const * const sum_dd_hat,
	int    const * const N,
	int    const * const k,
	double const * const Sigma_dd,
	double       * const Theta,
	double       * const D,
	int    const * const p,
	double       * const Transformation,
	double const * const minV,
	int const * const dl, double * dp, int * ip	)
{pfda_debug_step
	pfda_debug_cdl(debugnum_m5_0){
		pfda_debug_msg("pfda_m5_1 - \n");
		pfda_debug_msg("N:\n%d\n\n",*N );
		pfda_debug_argmat(delta, *N, *k);
		pfda_debug_arg(*k);
		pfda_debug_argmat(Theta, *p, *k);
		fflush(stdout);
	}
	/* memory allocation  and initialization */
	int obsnum;
	int k2=*k**k;
	int p2=*p**p;


	double * TST       = pfdaAlloc_d(p2,&dp);
	{//braces issolate the use of Sigma_d
	double inverseN=1/((double)(*N));
	double * Sigma_d=Transformation; ///Sigma_d borrows the memory from Transofrmation to use compute hold the sum, and related computations.
	pfda_debug_cdl(debugnum_m5_0){	pfda_debug_msg("sum conditional expectations\n");fflush(stdout);}
	dcopy_(&k2,sum_dd_hat, &one, Sigma_d,&one);
	pfda_debug_cdl(debugnum_m5_0){
		pfda_debug_argmat(Sigma_d, *k, *k);
		pfda_debug_arg(inverseN);
		fflush(stdout);
	}
	dscal_(&k2, &inverseN, Sigma_d, &one);
	pfda_matrix_outer_quadratic_form(TST,Theta,p,p,Sigma_d,k,dl, dp);
	pfda_debug_cdl(debugnum_m5_0){
		pfda_debug_argmat(Sigma_d, *k, *k);
		pfda_debug_msg("Compute eigen values\n");
		pfda_debug_argmat(TST, *p, *p);
		fflush(stdout);
		}
	Sigma_d=NULL;
	}

	double * vectors   = pfdaAlloc_d(p2,&dp);
	double * values    = pfdaAlloc_d(*p,&dp);
	pfda_eigens(TST,p,vectors, values,k,dl, dp, ip );

	pfda_debug_cdl(debugnum_m5_0){	pfda_debug_msg("compute transformation\n");fflush(stdout);	}
	/* delta_i<-Q^T\Theta\delta_i */
	dgemm_(&Trans,&NoTrans, k, k, p, &dOne, vectors, p, Theta, p, &dzero, Transformation, k);
	pfda_debug_cdl(debugnum_m5_0){
		pfda_debug_argmat(Transformation, *k, *k);;fflush(stdout);
		pfda_debug_msg("transform delta\n");fflush(stdout);	}

	double * delta_tmp = pfdaAlloc_d(*k,&dp);
	for(obsnum=0;obsnum<*N;obsnum++){
		dcopy_(k,delta+obsnum,N,delta_tmp,&one);
		dgemv_(&NoTrans, k, k, &dOne, Transformation, k, delta_tmp, &one, &dzero, delta+obsnum, N);
	}
	pfda_debug_cdl(debugnum_m5_0){pfda_debug_msg("copy new values into Theta and D\n");fflush(stdout);}
	int kp=*k**p;
	dcopy_(&kp, vectors, &one, Theta, &one);
	dcopy_(k, values, &one, D, &one);
	setMinVar(D, k, minV);

	pfda_debug_cdl(debugnum_m5_0){pfda_debug_msg("leaving pfda_m5_1\n");fflush(stdout);	}
}

/*!  Wrapper function that addes the necessary sumation of \f$ \sum\limits_{i=1}^N\widehat{\delta_i\delta_i^T} \f$ for the linear case.
\sa pfda_m5_0()
*/
void pfda_m5_1(
	int const * const N,
	double * const delta,
	int const * const k,
	double const * const Sigma_dd,
	double * const Theta,
	double * const D,
	int const * const p,
	double * const Transformation,
	double const * const minV,
	int const * const dl, double * dp, int * ip	)
{pfda_debug_step
	double * const sum_dd_hat = pfdaAlloc_d(*k**k,&dp);
	pfda_sum_cond_dd(sum_dd_hat, N, delta, Sigma_dd, k, dl);
	pfda_m5_0( delta, sum_dd_hat, N, k, Sigma_dd, Theta, D, p, Transformation, minV, dl, dp, ip);
}

/* Step 5 - estimates D_alpha and orthogonalizes Theta_f

@MEMORY
	dp length = k^2 + k + 2*p^2 + p +  max( k*p , 8*p)
*/
void pfdaSingle_m5(
	int    const * const N,
	double       * const delta,
	int    const * const k,
	double const * const Sigma_dd,
	double       * const Theta,
	double       * const D,
	int    const * const p,
	double const * const minV,
	int const * const dl, double * dp,int * ip
	)
{pfda_debug_step
	pfda_debug_cdl(debugnum_m5){pfda_debug_msg("pfdaSingle_m5\n");fflush(stdout);}
	double * transformation = pfdaAlloc_d(*k**k,&dp);
	pfda_m5_1(N, delta, k, Sigma_dd, Theta, D, p, transformation, minV, dl, dp, ip);
	pfdaFree(transformation);
	pfda_debug_cdl(debugnum_m5){pfda_debug_msg("leaving pfdaSingle_m5\n");fflush(stdout);}
}

/*!  Performs the E step for the univariate version.
@MEMORY
	- dp length= M + M*ka + N*ka + 2* ka^2
	- ip length = ka
*/
void pfdaSingle_e(
	double       * const alpha,
	double       * const Saa,
	double const * const y,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	int    const * const ka,
	double const * const B,
	int    const * const p,
	double const * const tm,
	double const * const tf,
	double const * const Da,
	double const * const sigma,
	const int * const dl, double *dp , int * ip)
{pfda_debug_step
	pfda_debug_cdl(debugnum_single_e){pfda_debug_msg("pfdaSingle_e - \n");fflush(stdout);}

	///* memory allocation */
	double * Ry       = pfdaAlloc_d(*M,&dp);
	double * phi      = pfdaAlloc_d(*M**ka, &dp);
	//double * alphaTmp = pfdaAlloc_d(*ka**N,&dp);
	double inv_epsilon = 1.0 / *sigma;
	pfda_debug_cdl(debugnum_single_e){pfda_debug_msg("sigma epsilon inverse: %g\n",inv_epsilon);fflush(stdout);}

	pfda_computeResid(Ry, y, nobs, M, N, ka, B, p, tm, NULL, NULL, dl, dp);
	// the nulls mean that the PC functions are not subtracted off and the phi matrix is not computed.   The memoy requirement for the above is 0.
	pfda_debug_cdl(debugnum_single_e){pfda_debug_msg("Ry:\n");
				printyvec(Ry, nobs,N);fflush(stdout);}

	dgemm_(&NoTrans, &NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);//compute BTf = B*tf
	pfda_debug_cdl(debugnum_single_e){pfda_debug_msg("phi:\n");
				printmat(phi, *M, *ka);fflush(stdout);}
	double *phii=phi, *Sigma_i=Saa, *Ryi=Ry;
	for(int obsnum=0;obsnum<*N;obsnum++){
		pfdaSingle_e_1( alpha+obsnum, N, Sigma_i, Ryi, nobs+obsnum, M, ka, phii, Da, &inv_epsilon, dl, dp, ip);
		phii+=nobs[obsnum];
		Ryi +=nobs[obsnum];
		Sigma_i += *ka**ka;
	}
}

/*!  Step for E step
	The inner workings of the pfdaSingle_e step but for only a single subject at a time.
MEMORY:
	dp length = ka^2 + ka
	ip length = ka
*/
void pfdaSingle_e_1(
	double       * const alpha_i,
	int    const * const ldalpha,
	double       * const Sigma_i,
	double const * const Ryi,
	int    const * const nobs_i,
	int    const * const M,
	int    const * const ka,
	double const * const phii,
	double const * const Da,
	double const * const inv_epsilon,
	int const * const dl, double *dp , int * ip)
{pfda_debug_step
	double * alphaTmp = pfdaAlloc_d(*ka,&dp);
	pfda_debug_cdl(debugnum_single_e_inloop){pfda_debug_argmat2(*nobs_i,*ka,phii,M);fflush(stdout);}

	dsyrk_(&Upper, &Trans, ka, nobs_i, inv_epsilon, phii, M, &dzero, Sigma_i, ka);
	pfda_debug_cdl(debugnum_single_e_inloop){pfda_debug_msg("phii^t phii:\n");pfda_debug_argmat(Sigma_i,*ka,*ka);fflush(stdout);}
	addInvDiag(ka, Da, Sigma_i, ka);
	pfda_debug_cdl(debugnum_single_e_inloop){pfda_debug_msg("Da^-1+phi^t phi:\n");pfda_debug_argmat(Sigma_i,*ka,*ka);fflush(stdout);}

	int sr=0;
	pfda_sym_inverse(Sigma_i,ka, &sr, dl, dp, ip);
	if(sr>0) pfda_error("PFDA ERR - pfdaSingle_e: pfda_sym_inverse returned %d\n",sr);
	pfda_debug_cdl(debugnum_single_e_inloop){
		pfda_debug_msg("pfda_sym_inverse Sigma_i success\n");
		pfda_debug_argmat(Sigma_i,*ka,*ka);fflush(stdout);
	}

	dgemv_(&Trans, nobs_i, ka, &dOne, phii, M, Ryi, &one, &dzero, alphaTmp, &one); 	// alphaTmp = phii^T * Ry
	pfda_debug_cdl(debugnum_single_e_inloop){pfda_debug_msg("phii^t Ryi:\n");pfda_debug_argmat2(one,*ka,alphaTmp,&one);fflush(stdout);}
	dsymv_(&Upper, ka, inv_epsilon, Sigma_i, ka, alphaTmp, &one, &dzero, alpha_i, ldalpha);		// alpha = 1/sigma * Sigma_i * alphaTmp
	pfda_debug_cdl(debugnum_single_e_inloop){pfda_debug_argmat2(one,*ka,alpha_i,ldalpha);fflush(stdout);}
}

/*! performs the univariate fit for original model of functional data analysis using functional principle components
INPUTS:
	t		Mx1 vector of time values
	y		Mx1 vector of y values   (Also an output for residuals)
	nobs		Nx1 vector of number of observations per subject
	M		number of observations
	N		number of subjects
	ka		number of principle components used for modeling y (ka <= p)
	B		Mxp matrix of timepoints evaluated at the orthogonalized basis functions
	p		number of orthogonal basis functions
	minV		for ridge regression type computations
	--not in find inits--
	penaties	4 vector of penaltie values in the form 	lambda_mu, lambda_nu, lambda_f, lambda_g
	K		pxp matrix of the integrated square second derivative of the basis functions
OUTPUTS:	name|Description(name in R)
	y		residuals
	tm		px1 vector estimating mean structore for y (theta_mu)
	tf		pxka matrix representing the subject specific portion of the estimate for y(Theta_f)
	Da		ka vector representing the variances of alpha (D_alpha)
	sigma		variance of y residuals (sigma_epsilon)
	alpha		Nxka matrix of PC coefficients of y (alpha)
	Saa		kaxkaxN 3 dimentional array of variance components (Sigma_alpha_alpha)
CONTROL VALUES:
	tol		Tolerance for determining convergence
	maxI	on Input: the Maximum number of iterations.
			on output: THe numer of itterations used.
	dl	controls printing of debug information
	incInits	Boolean: does the passing in information include initial values? if true findInits is not called
@MEMORY
	dp	length = N*p^2 + p + p*ka + ka +
		MAX:
			1. step 1                  	= M + M*ka + 2*ka^2   (exceeded by E)
			2. step 2                 	= p^2+ M+ M*ka
			3. step 3                 	= M + p^2 + p + ka^2
			5. step 4/5             	= ka+ 2*p^2 + p +  max(  ka, 8 ) * p
			E. pfdaSingle_e       	= M + M*ka + N*ka + 2* ka^2
			others              	= M
	ip	length = 8 p
*/
void pfdaSingle(
	/* input Values */
	double const * const t,
	double * const y, 		//on exit y will hold residuals.
	//int const * const  obs,
	int const * const  nobs,
	int const * const  M,
	int const * const  N,
	int const * const  ka,
	double const * const B,
	int const * const  p,
	double const * const minV,
	double const * const lm,
	double const * const lf,
	double const * const K,
	/* State Values */
	double * const tm,
	double * const tf,
	double * const Da,
	double * const sigma,
	double * const alpha,
	double * const Saa,
	/* Control Values */
	double * const tol,
	int * const maxI,
	int * const dl,
	double *dp , int * ip)
{pfda_debug_step
if(*dl){
pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Entering pfdaSingle\n");fflush(stdout);}
pfda_debug_cdl(debugnum_single_inputs){
	pfda_debug_msg("INPUTS\n----------------------------------------\n");
	pfda_debug_arg(*M);
	pfda_debug_msg("N:\t%d\n",*N );
	pfda_debug_arg(*ka);
	pfda_debug_arg(*p);
	pfda_debug_arg(*lm);
	pfda_debug_arg(*lf);
	pfda_debug_msg("sigma:\n%g\n", *sigma);
	pfda_debug_arg(*minV);
	pfda_debug_msg("Tolerance:\t%g\nMaximun Iterations:\t%d\n",  *tol, *maxI);
	fflush(stdout);

	pfda_debug_msg("Vector Arguments:\n");
	pfda_debug_argmat(t, 1, *M);
	pfda_debug_argmat(y, 1, *M);
	pfda_debug_msg("nobs:\n"); for(int i=0;i<*N;i++)pfda_debug_msg("%d\t",nobs[i]); pfda_debug_msg("\n\n");
	fflush(stdout);

	pfda_debug_msg("Outputs:");
	pfda_debug_argmat(tm, 1, *p);
	pfda_debug_argmat(tf, *p, *ka);
	pfda_debug_argmat(Da, 1, *ka);
	pfda_debug_argmat(alpha, *N, *ka);
	fflush(stdout);
}
pfda_debug_cdl(debugnum_single_inputs_large){
	pfda_debug_argmat(K, *p, *p);
	pfda_debug_argmat(B, *M, *p);
	pfda_debug_argmat(Saa, *N, *ka**ka);
}
}
pfdaFindInitSingle(t, y, nobs, M, N, ka,  B,  p,  minV, //input values
			tm,  tf,  Da,  sigma, alpha, dl); //return values
setMinVar(Da, ka, minV);
setMinVar(sigma, &one, minV);

double * btb = pfdaAlloc_d(*p**p**N,&dp);
pfda_computebtb(btb,N,B,M,p,nobs,dl);

int I=0;
double   sigma_old=0;
double   convergenceCriteria=0;
double * tmOld  = pfdaAlloc_d(*p    , &dp);
double * tf_old = pfdaAlloc_d(*p**ka, &dp);
double * Da_old = pfdaAlloc_d(*ka   , &dp);
while(I < *maxI)/*  */{
{///* setup for convergence */
	pfda_debug_cdl(debugnum_singe_steps){
		pfda_debug_msg("Entering loop: %d\n",I);
		pfda_debug_msg("Setup for convergence\n");
		fflush(stdout);
	}
	sigma_old = *sigma; //teporarily store values for use with convergence
	dcopy_(p, tm, &one, tmOld, &one);
	int pka=*p**ka;dcopy_(&pka, tf, &one, tf_old, &one);
	dcopy_(ka, Da, &one, Da_old, &one);
}
{///* step 1 */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Step 1\n");fflush(stdout);}
	pfda_m1(sigma,y,nobs,M,N,ka,B,p,minV,tm,tf,alpha,Saa,dl, dp);
	pfda_debug_cdl(debugnum_single_estimates){pfda_debug_msg("sigma: %# .8g\n\n", *sigma);fflush(stdout);}
}
{///* step 2 */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Step 2\n");fflush(stdout);}
	pfda_m2(tm,y,nobs,M,N,ka,B,p,lm,K,tf,sigma,alpha,dl, dp);
	pfda_debug_cdl(debugnum_single_estimates){pfda_debug_argmat(tm,*p,one);fflush(stdout);}
}
{///* step 3 */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Step 3\n");fflush(stdout);}
	pfda_m3(tf, y, nobs, M, N, ka, B, p, lf, K, tm, sigma, alpha, Saa, btb, dl, dp, ip);
	pfda_debug_cdl(debugnum_single_estimates){pfda_debug_argmat(tf,*p,*ka);fflush(stdout);}
}
{///* Steps 4 & 5 */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Step 4&5\n");fflush(stdout);}
	pfdaSingle_m5(N, alpha, ka, Saa, tf, Da, p, minV, dl, dp, ip);
	pfda_debug_cdl(debugnum_single_estimates){
		pfda_debug_argmat(tf,*p,*ka);
		pfda_debug_argmat(Da,*ka,one);
		pfda_debug_argmat(alpha,*p,*ka);
		fflush(stdout);}
}
{///* E-Step */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Step E\n");fflush(stdout);}
	pfdaSingle_e(alpha, Saa, y, nobs, M, N, ka, B, p, tm, tf, Da, sigma, dl, dp, ip);
	pfda_debug_cdl(debugnum_single_estimates){
		pfda_debug_argmat(alpha,*p,*ka);
		pfda_debug_argmat(Saa,*N,*ka**ka);
		fflush(stdout);}
}
{///* Compute Convergence Criteria */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Computing convergence Criteria\n");fflush(stdout);}
	convergenceCriteria = fabs((sigma_old-*sigma)/(*sigma));

	daxpy_(p, &mOne, tm, &one, tmOld, &one);
	double cctm=0;
	for(int i=0;i<*p;i++){
		cctm+=fabs(tmOld[i]/tm[i]);
	}
	convergenceCriteria+=cctm;

	double cctf=0;
	{int pka=*p**ka;daxpy_(&pka, &mOne, tf, &one, tf_old, &one);
	for(int i=0;i<pka;i++)cctf+=fabs(tf_old[i]/tf[i]);
	convergenceCriteria+=cctf;}

	daxpy_(ka, &mOne, Da, &one, Da_old, &one);
	double ccDa=0;
	for(int i=0;i<*ka;i++)ccDa+=convergenceCriteria+=fabs(Da_old[i]/Da[i]);
	convergenceCriteria+=ccDa;

	///* Debug */
	pfda_debug_cdl(debugnum_single_convergence_components){
			pfda_debug_msg("Convergence criteria components:\n");
			pfda_debug_msg("sigma: \t%5.5g\n", fabs(sigma_old-*sigma));
			pfda_debug_msg("tm:   \t%5.5g\n", cctm);
			pfda_debug_msg("Tf:   \t%5.5g\n", cctf);
			pfda_debug_msg("Da:   \t%5.5g\n", ccDa);
			fflush(stdout);
	}
	pfda_debug_cdl(debugnum_single_convergence){pfda_debug_msg("Criteria:%g\n", convergenceCriteria);fflush(stdout);}
}
I++;
if(convergenceCriteria < *tol)break;
}//end while loop
{	///* Finishing */
pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Finishing\n");fflush(stdout);}
pfda_computeResid(y, y, nobs, M, N, ka, B, p, tm, tf, alpha, dl, dp);
*tol = convergenceCriteria;
*maxI = I;
}
pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Leaving pfdaSingle\n");fflush(stdout);}
}//end pfdaSingle


//---------------------------------------------------
// this marks the new section to allow for the extension of
// unpenalized variables included in the model


/*! Estimates alpha, Saa, and aa
 *
 */
void single_c_E(
	double * const alpha,
	double * const aa,
	double * const Saa,
	const double * const y,
	const double * const Z,
	const double * const B,
	const double * const tz,
	const double * const tm,
	const double * const tf,
	const double * const Da,
	const double * const sigma,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const k,
	const int * const p,
	const int * const dl,
	double * dp, int*ip)
{
	double * Ry = pfdaAlloc_d(*M,&dp);
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(Z && kz && tz && *kz)dgemv_(&NoTrans, M, kz, &mOne, Z, M, tz, &one, &dOne, Ry, &one); //y = (y-Z*tz)
	pfdaSingle_e(alpha, Saa, Ry, nobs, M, N, k, B, p, tm, tf, Da, sigma, dl, dp, ip);
}

/*! Fits the unpenalized portion of the model
 *
 *	@MEMORY
 *	-dp = M + kz^2
 */
void single_c_unpenalized(
	double       * const tz,
	const double * const y,
	const double * const Z,
	const double * const B,
	const double * const tm,
	const double * const tf,
	const double * const alpha,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const k,
	const int * const p,
	const int * const dl,
	double * dp, int*ip)
{
	double * Ry = pfdaAlloc_d(*M, &dp);
	single_c_resid( Ry, y, NULL, nobs, M, N, NULL, k, B, p, NULL, tm, tf, alpha, dl, dp);

	if(Z && *kz){
		double * left = pfdaAlloc_d(*kz**kz,&dp);
		dsyrk_(&Upper, &Trans, kz, M, &dOne, Z, M, &dzero, left, kz);
		dgemv_(&Trans, M, kz, &dOne, Z, M, Ry, &one, &dzero, tz, &one);
		int sr =0, lwork = 10**kz;
		dsysv_(&Upper, kz, &one, left, kz, ip, tz, kz, dp, &lwork, &sr);
		if(sr)pfda_error("PFDA ERROR:dual_ca_i: computing tz dsysv returned:%d",sr);
	}
}

/*! Fits penalized mean function
 *
 * @MEMORY
 * -dp
 * -lp
 */
void single_c_penalized(
	double       * const tm,
	const double * const y,
	const double * const Z,
	const double * const B,
	const double * const tz,
	const double * const tf,
	const double * const alpha,
	const double * const sigma,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const k,
	const int * const p,
	const double * const lm,
	const double * const K,
	const int * const dl,
	double * dp, int*ip)
{
	double * Ry = pfdaAlloc_d(*M,&dp);
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(Z && kz && tz && *kz)dgemv_(&NoTrans, M, kz, &mOne, Z, M, tz, &one, &dOne, Ry, &one); //y = (y-Z*tz)
	pfda_m2(tm, Ry, nobs, M, N, k, B, p, lm, K, tf, sigma, alpha, dl, dp);
}

/*! Fits the principle components portion of theunivariate model
 *
 */
void single_c_princcomp(
	double * const tf,
	const double * const y,
	const double * const Z,
	const double * const B,
	const double * const tz,
	const double * const tm,
	const double * const alpha,
	const double * const sigma,
	const double * const aa,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const k,
	const int * const p,
	const double * const lf,
	const double * const K,
	const double * const btb,
	const int * const dl,
	double * dp, int*ip)
{
	double * Ry = pfdaAlloc_d(*M,&dp);
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(Z && kz && tz && *kz)dgemv_(&NoTrans, M, kz, &mOne, Z, M, tz, &one, &dOne, Ry, &one); //y = (y-Z*tz)
	// double * btb = pfdaAlloc_d(*p**p**N,&dp);pfda_computebtb(btb,N,B,M,p,nobs,dl);
	pfda_m3_core( tf, aa, Ry, nobs, M, N, k, B, p, lf, K, tm, sigma, alpha, btb, dl, dp, ip);
}

/*! estimates variances for the single model
 *
 */
void single_c_variances(
	      double * const tf,
	      double * const alpha,
	      double * const Da,
	      double * const sigma,
	const double * const y,
	const double * const Z,
	const double * const B,
	const double * const tz,
	const double * const tm,
	const double * const Saa,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const k,
	const int * const p,
	const int * const dl,
	double * dp, int* ip)
{
	double * Ry = pfdaAlloc_d(*M,&dp);
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(Z && kz && tz && *kz)dgemv_(&NoTrans, M, kz, &mOne, Z, M, tz, &one, &dOne, Ry, &one); //y = (y-Z*tz)
	pfda_m1(sigma,y,nobs,M,N,k,B,p,&dzero,tm,tf,alpha,Saa,dl, dp);
	pfdaSingle_m5(N, alpha, k, Saa, tf, Da, p, &dzero, dl, dp, ip);
}

/*! complete single algorithm that also factors in Z
 *
 */
void single_c_core(
	double * const y,
	double * const Z,
	double * const B,
	double * const tz,
	double * const tm,
	double * const tf,
	double * const alpha,
	double * const Da,
	double * const sigma,
	double * const aa,
	double * const Saa,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const k,
	const int * const p,
	const double * const lm,
	const double * const lf,
	const double * const K,
	const double * const minV,
	      int * const maxI,
	      double * const tol,
	const int * const dl, double * dp, int * ip)
{
if(*dl){
pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Entering pfdaSingle\n");fflush(stdout);}
pfda_debug_cdl(debugnum_single_inputs){
	pfda_debug_msg("INPUTS\n----------------------------------------\n");
	pfda_debug_arg(*M);
	pfda_debug_arg(*N);
	pfda_debug_arg(*k);
	pfda_debug_arg(*p);
	pfda_debug_arg(*lm);
	pfda_debug_arg(*lf);
	pfda_debug_arg(*sigma);
	pfda_debug_arg(*minV);
	pfda_debug_msg("Tolerance:\t%g\nMaximun Iterations:\t%d\n",  *tol, *maxI);
	fflush(stdout);

	pfda_debug_msg("Vector Arguments:\n");
	pfda_debug_yvec(y);
	pfda_debug_argveci(nobs,N);
	fflush(stdout);

	pfda_debug_msg("Outputs:");
	pfda_debug_argmat(tm, 1, *p);
	pfda_debug_argmat(tf, *p, *k);
	pfda_debug_argmat(Da, 1, *k);
	pfda_debug_argmat(alpha, *N, *k);
	fflush(stdout);
}
pfda_debug_cdl(debugnum_single_inputs_large){
	pfda_debug_argmat(K, *p, *p);
	pfda_debug_argmat(B, *M, *p);
	pfda_debug_argmat(Saa, *N, *k**k);
}
}

double * btb = pfdaAlloc_d(*p**p**N,&dp);
pfda_computebtb(btb,N,B,M,p,nobs,dl);
pfda_s_i(tm, tf, alpha, Da, aa, sigma, y, nobs, M, N, k, B, btb, p, minV, dl, dp, ip);
int I=0;
double   sigma_old=0;
double   convergenceCriteria=0;
double * tmOld  = pfdaAlloc_d(*p    , &dp);
double * tf_old = pfdaAlloc_d(*p**k, &dp);
double * Da_old = pfdaAlloc_d(*k   , &dp);
while(I < *maxI){
	{///* setup for convergence */
		pfda_debug_cdl(debugnum_singe_steps){
			pfda_debug_msg("Entering loop: %d\n",I);
			pfda_debug_msg("Setup for convergence\n");
			fflush(stdout);
		}
		sigma_old = *sigma; //teporarily store values for use with convergence
		dcopy_(p, tm, &one, tmOld, &one);
		int pka=*p**k;dcopy_(&pka, tf, &one, tf_old, &one);
		dcopy_(k, Da, &one, Da_old, &one);
	}
	single_c_E( alpha, aa, Saa, y, Z, B, tz, tm, tf, Da, sigma, nobs, N, M, kz, k, p, dl, dp, ip);
	single_c_unpenalized( tz, y, Z, B, tm, tf, alpha, nobs, N, M, kz, k, p, dl, dp, ip);
	single_c_penalized( tm, y, Z, B, tz, tf, alpha, sigma, nobs, N, M, kz, k, p, lm, K, dl, dp, ip);
	single_c_princcomp( tf, y, Z, B, tz, tm, alpha, sigma, aa, nobs, N, M, kz, k, p, lf, K, btb, dl, dp, ip);
	single_c_variances( tf, alpha, Da, sigma, y, Z, B, tz, tm, Saa, nobs, N, M, kz, k, p, dl, dp, ip);
	for(int i=0;i<*k;i++){if(Da[i]<*minV)Da[i]=*minV;} /// Enforces Minimun variances on \f$D_\alpha\f$
	{// Compute Convergence Criteria
		pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Computing convergence Criteria\n");fflush(stdout);}
		convergenceCriteria = fabs((sigma_old-*sigma)/(*sigma));

		daxpy_(p, &mOne, tm, &one, tmOld, &one);
		double cctm=0;
		for(int i=0;i<*p;i++){
			cctm+=fabs(tmOld[i]/tm[i]);
		}
		convergenceCriteria+=cctm;

		double cctf=0;
		{int pka=*p**k;daxpy_(&pka, &mOne, tf, &one, tf_old, &one);
		for(int i=0;i<pka;i++)cctf+=fabs(tf_old[i]/tf[i]);
		convergenceCriteria+=cctf;}

		daxpy_(k, &mOne, Da, &one, Da_old, &one);
		double ccDa=0;
		for(int i=0;i<*k;i++)ccDa+=convergenceCriteria+=fabs(Da_old[i]/Da[i]);
		convergenceCriteria+=ccDa;

		pfda_debug_cdl(debugnum_single_convergence_components){
				pfda_debug_msg("Convergence criteria components:\n");
				pfda_debug_msg("sigma: \t%5.5g\n", fabs(sigma_old-*sigma));
				pfda_debug_msg("tm:   \t%5.5g\n", cctm);
				pfda_debug_msg("Tf:   \t%5.5g\n", cctf);
				pfda_debug_msg("Da:   \t%5.5g\n", ccDa);
				fflush(stdout);
		}
		pfda_debug_cdl(debugnum_single_convergence){pfda_debug_msg("Criteria:%g\n", convergenceCriteria);fflush(stdout);}
	}
	I++;
	if(convergenceCriteria < *tol)break;
}//end while loop
{	///* Finishing */
	pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Finishing\n");fflush(stdout);}
	if(*maxI<=I){pfda_warning("EM-algorithm did not converge");}
	single_c_resid( y, y, Z, nobs, M, N, kz, k, B, p, tz, tm, tf, alpha, dl, dp);
	*tol = convergenceCriteria;
	*maxI = I;
}
pfda_debug_cdl(debugnum_singe_steps){pfda_debug_msg("Leaving single_c_core\n");fflush(stdout);}
}//end pfdaSingle








