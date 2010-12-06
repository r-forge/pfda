/*! bin_E(W).c
coputes the expected values for E(W|Y,Z) and E(W|Y) by simulation methods
*/
#include "pfda_error.h"
#include "pfda_gen.h"
#include "pfda_memory.h"
#include "pfda_single.h"
#include "pfda_utils.h"
#include <R_ext/BLAS.h>

/*! 	takes a vector old and copies to new without row droprow. 
	works for vectors by specifying ncol = 1
	
	\ingroup utils
*/
void pfda_matrix_drop1(
	double       * const new,
	int    const * const ldnew,
	double const * const old,
	int    const * const ldold,
	int    const * const nrow,
	int    const * const ncol,
	int    const * const droprow,
	int const * const dl)
{pfda_debug_step
	// pfda_debug_cdl(debugnum_pfda_matrix_drop1)pfda_debug_msg("Entering pfda_matrix_drop1\n");
	int col=0;
	for(col=0;col<*ncol;col++){
		pfda_debug_line;
		int newrow=0, oldrow=0;
		while(oldrow<*nrow)if(oldrow!=*droprow)new[col**ldnew + newrow++]=old[col**ldold + oldrow++]; else oldrow++;
	}
	// pfda_debug_cdl(debugnum_pfda_matrix_drop1)pfda_debug_msg("exiting pfda_matrix_drop1\n");
}

/*!  computes the parameters for generating W|Y
\ingroup binary_single
@MEMORY
	- dp length = 3 k 
	- ip length = k
*/
void pfda_bin_single_generate_w_parms1(
	double * const mu,
	double * const sigma,
	double const * const RW_trimmed,
	double const * const Da,
	double const * const phi_trimmed,
	int    const * const k,
	int    const * const ni_trimmed,
	int const * const dl,
	double * dp, int * ip)
{pfda_debug_step
	// pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms1)pfda_debug_msg("Entering debugnum_pfda_bin_single_generate_w_parms1\n");
	// /* compute \mu^* and \sigma^* */
	pfdaSingle_e_1( mu, &one, sigma, RW_trimmed, ni_trimmed, ni_trimmed, k, phi_trimmed, Da, &dOne, dl, dp, ip);
	// pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms1)pfda_debug_msg("Exiting debugnum_pfda_bin_single_generate_w_parms1\n");
}

/*! 	Computes parameters for generating W|Y

\ingroup binary_single
@MEMORY
	- dp length = p + p^2 + k*p
*/
void pfda_bin_single_generate_w_parms2(
	double       * const a,
	double       * const s,
	double const * const mu,
	double const * const sigma,
	double const * const B,
	int    const * const M, 
	int    const * const p, 
	double const * const tm,
	double const * const tf,
	int    const * const k,
	int const * const dl,
	double * dp)
{pfda_debug_step
	pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms2){
		pfda_debug_msg("Entering pfda_bin_single_generate_w_parms2\n");
		pfda_debug_msg("tf:\n");printmat(tf,*p,*k);
		pfda_debug_msg("mu:\n");printmat(tf,*k,one);
		
		fflush(stdout);
	};
	// /* a_ij */
	pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms2)pfda_debug_msg("Computing a\n");
	double * tfmu = pfdaAlloc_d(*p, &dp);
	dgemv_(&NoTrans, p,k, &dOne, tf, p, mu, &one, &dzero, tfmu, &one);
	daxpy_(p, &dOne, tm, &one, tfmu, &one);
	*a = ddot_(p, B, M, tfmu, &one);
	pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms2){pfda_debug_msg("a:  %g\n", *a);fflush(stdout);}
	
	// /* sigma_ij */
	double * tfb = pfdaAlloc_d(*k,&dp);
	dgemv_(&Trans, p, k, &dOne, tf, p, B, M, &dzero, tfb,&one);
	pfda_matrix_inner_quadratic_form(s,tfb,&one,k,sigma,k,dl,dp);
	
	
	// double * tmp = pfdaAlloc_d(*p**p,&dp);
	// pfda_matrix_outer_quadratic_form(tmp,tf,p,p,sigma,k,dl,dp );
	// pfda_matrix_inner_quadratic_form(s,B,&one,M,tmp,p,dl,dp );
	*s += dOne;
	pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms2){pfda_debug_msg("s:  %g\n", *s);fflush(stdout);}
	pfda_debug_cdl(debugnum_pfda_bin_single_generate_w_parms2)pfda_debug_msg("Exiting pfda_bin_single_generate_w_parms2\n");
}


/*! 	Generates a set of simulationed version for a single observation of W.

\ingroup binary_single
@MEMORY
	- dp length = ni_trimmed * k + ni_trimmed + k^2 + k + p + p^2 + k*p + 3 k 
	- ip length = k
*/
void pfda_bin_s_gen_w(
	double       * const w_sim,
	int    const * const Nsim,
	int    const * const Yi, 
	double const * const RWi,
	int    const * const ni,
	double const * const Bi,
	int    const * const M, 
	int    const * const p, 
	double const * const tm,
	double const * const tf,
	int    const * const k,
	double const * const Da,
	int    const * const j,
	int const * const dl,
	double * dp, int * ip)
{
	pfda_debug_cdl(debugnum_pfda_bin_s_gen_w)pfda_debug_msg("entering pfda_bin_s_gen_w\n");
	pfda_debug_cdl(debugnum_pfda_bin_s_gen_w){
		pfda_debug_msg("RWi: \n");printmat(RWi,*ni,one);
		fflush(stdout);
	}	
	int ni_trimmed = (*ni-1);
	double * phi_trimmed = pfdaAlloc_d(ni_trimmed**k, &dp);
	double * RW_trimmed  = pfdaAlloc_d(ni_trimmed, &dp);
	
	// /* setup trimmed variables*/
	for(int l=0;l<*ni;l++){
		if(l!=*j){
			int o=l;
			if(l>*j)o-=1;
			dgemv_(&Trans, p,k, &dOne, tf, p, Bi+l, M, &dzero, phi_trimmed+o, &ni_trimmed);
			RW_trimmed[o] = RWi[l]; 
		}
	}
	pfda_debug_cdl(debugnum_pfda_bin_s_gen_w){
		pfda_debug_msg("RW_trimmed: \n");printmat(RW_trimmed,ni_trimmed,one);
		fflush(stdout);
	}
	
	// /* compute parameters */
	double * sigma = pfdaAlloc_d(*k**k, &dp);
	double * mu    = pfdaAlloc_d(*k, &dp);
	double a_ij=0.0, sigma_ij=0.0;
	pfda_bin_single_generate_w_parms1( mu, sigma, RW_trimmed, Da, phi_trimmed, k, &ni_trimmed, dl, dp, ip);
	pfda_bin_single_generate_w_parms2( &a_ij, &sigma_ij, mu, sigma, Bi+*j, M, p, tm, tf, k, dl, dp);
	pfda_debug_cdl(debugnum_pfda_bin_s_gen_w){
		pfda_debug_msg("mu: \n");printmat(mu,*k,one);
		pfda_debug_msg("sigma: \n");printmat(sigma,*k,*k);
		pfda_debug_msg("a_ij: %g\n",a_ij);
		pfda_debug_msg("sigma_ij: %g\n",sigma_ij);
		fflush(stdout);
	}
	// /* Generate Simulated Variables */
	if(Yi[*j]){
		double c = -a_ij/sigma_ij;
		pfda_gen_truncnorm( w_sim, Nsim, &c, dl);
		for(int k=0;k<*Nsim;k++){
			w_sim[k] *= sigma_ij;
			w_sim[k] += a_ij;
		}
		// pfda_rtruncnorm(w_sim, Nsim, &dzero, &one, &dzero, NULL, &a_ij,	&one, &sigma_ij,	&one, dl);
	} else {
		double c = a_ij/sigma_ij;
		pfda_gen_truncnorm( w_sim, Nsim, &c, dl);
		for(int k=0;k<*Nsim;k++){
			w_sim[k] *= -sigma_ij;
			w_sim[k] += a_ij;
		}
		// pfda_rtruncnorm(w_sim, Nsim, &dzero, NULL, &dzero, &one, &a_ij,	&one, &sigma_ij,	&one, dl);
	}
}

/*!	interior of for loop separated into separate function for memory management reasons.

\ingroup binary_single
@MEMORY:
	- dp length 	= ni*kr + kr + (ni_trimmed * k + ni_trimmed + k^2 + k + p + p^2 + k*p + 3 k)
	- ip length = k
*/
void pfda_bin_single_approximate_moments_forobs(
	double * const w_i,
	double * const ww_i,
	int    const * const y_i,
	int    const * const n_i,
	int    const * const M,
	int    const * const p,
	int    const * const k,
	double const * const tm,
	double const * const Bi,
	double const * const tf,
	double const * const Da,
	int    const * const kr,
	double const * const weight,
	int const * const dl, double * dp, int * ip)
{
	/// Rversion is almost .bin_s_w_1 but includes the weight.
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments_forobs){
		pfda_debug_msg("Entering pfda_bin_single_approximate_moments_forobs\n");
		pfda_debug_msg("Weight: %g\n",*weight);
	}
	// /* setup for loop */
	double * RW_i = pfdaAlloc_d(*n_i,&dp);
	pfda_computeResid(RW_i,w_i,n_i, M, NULL, k, Bi, p, tm, NULL, NULL, dl, dp);
	
	double * wsim = pfdaAlloc_d(*n_i**kr, &dp);
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments_forobs){
		pfda_debug_msg("wsim:\n");printmat(wsim,*kr,*n_i);
	}
	for(int j=0;j<*n_i;j++)
		pfda_bin_s_gen_w( wsim+*kr*j, kr, y_i, RW_i, n_i, Bi, M, p, tm, tf, k, Da, &j, dl, dp, ip);
		
	double krinv = 1.0/(double)*kr;
	double * V1 = pfdaAlloc_d(*kr,&dp);
	for(int i=0;i<*kr;i++)V1[i]=krinv;
	
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments_forobs){
		pfda_debug_msg("wsim:\n");printmat(wsim,*kr,*n_i);
	}
	double cweight = (1-*weight);
	dgemv_(&Trans, kr, n_i, weight, wsim, kr, V1, &one, &cweight, w_i, &one);

	double w2 = *weight/(double)*kr;
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments_forobs){
		pfda_debug_msg("w2: %g\n",w2);
		pfda_debug_msg("wsim:\n");printmat(wsim,*kr,*n_i);
	}
	dsyrk_(&Upper, &Trans, n_i, kr, &w2, wsim, kr, &cweight, ww_i, n_i);
	pfda_fillsym (ww_i,n_i,dl);
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments_forobs)
		pfda_debug_msg("Exiting pfda_bin_single_approximate_moments_forobs\n");
}

/*!	the final step in the Binary Paired fucntional data anlalysis paper.

	finds \f$ E(W|Y) \f$ and \f$ E(WW|Y) \f$ 
	
\ingroup binary_single
@MEMORY
	upper bounds
	- dp length = ni*kr + kr + (ni * k + ni + k^2 + k + p + p^2 + k*p + 3 k)
	where ni=max(nobs)
	- ip length = k
*/
void pfda_bin_single_approximate_moments(
	double       * const w,
	double       * const ww,
	int    const * const y,
	int    const * const nobs,
	int    const * const M,
	int    const * const N,
	double const * const B,
	int    const * const p,
	int    const * const k,
	double       * const tm,
	double       * const tf,
	double       * const Da,
	double const * const weight,
	int    const * const kr,
	int const * const dl,
	double * dp, int * ip)
{
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments){
		pfda_debug_msg("Entering pfda_bin_single_approximate_moments\n");
		pfda_debug_msg("w:\n");printyvec(w,nobs,N);
		pfda_debug_msg("ww:\n");double*wwi=ww;for(int i=0;i<*N;i++){pfda_debug_msg("ww[%d]:\n",i);printmat(wwi,nobs[i],nobs[i]);wwi+=nobs[i]*nobs[i];}
		pfda_debug_msg("y:\n");printyveci(y,nobs,N);
		pfda_debug_msg("nobs:\n");printmati(nobs,one,*N);
		pfda_debug_msg("M: %d\n\n",*M);
		pfda_debug_msg("N: %d\n\n",*N);
		pfda_debug_msg("p: %d\n\n",*p);
		pfda_debug_msg("k: %d\n\n",*k);
		pfda_debug_msg("tm:\n");printmat(tm,*p,one);
		pfda_debug_msg("tf:\n");printmat(tf,*p,*k);
		pfda_debug_msg("Da:\n");printmat(Da,one,*k);
		pfda_debug_msg("weight: %g\n\n",*weight);
		pfda_debug_msg("kr: %d\n\n",*kr);
		fflush(stdout);
	}
	int const * y_i = y;
	double * w_i=w, *ww_i = ww;
	double const *Bi=B;
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments){pfda_debug_msg("Entering obsnum loop\n");fflush(stdout);}
	for(int obsnum=0;obsnum<*N;obsnum++){
		pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments){pfda_debug_msg("obsnum= %d\n",obsnum);fflush(stdout);}
		pfda_bin_single_approximate_moments_forobs(w_i,ww_i, y_i, nobs+obsnum, M, p, k, tm, Bi, tf, Da, kr, weight, dl, dp, ip);
		w_i     += nobs[obsnum];
		ww_i    += nobs[obsnum]*nobs[obsnum];
		y_i     += nobs[obsnum];
		Bi   += nobs[obsnum];
	}
	pfda_debug_cdl(debugnum_pfda_bin_single_approximate_moments)pfda_debug_msg("leaving pfda_bin_single_approximate_moments\n");
}

