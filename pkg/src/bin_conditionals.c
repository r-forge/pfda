/*! \file bin_conditionals.c

by Andrew Redd
Performs the conditional expectation for the bivariate model with a binary varible.
*/
#include "pfda.h"

void pfda_bin_cond_aa(
	/* output  */
	double * const aa_hat,
	/* input */
	double const * const w_i, 
	double const * const ww_i, 
	int    const * const n,
	int    const * const M,
	double const * const S1_i,
	double const * const phi_i,
	double const * const Sigma_aa_i,
	int    const * const ka,
	
	/* Control Values */
	int const * const dl,
	double * dp)
/* DESCRIPTION
computes the  conditional expectation of \widehat{\alpha\alpha\trans}
INPUTS:
	w = the expected value of w given y and z.
	ww_i = the expected value of ww^t given y and z
	n = the number of observations for this subject/curve
	M = the total number of observations, number of rows of phi_i.
	S1_i = vector of length ka from the computations from pfdaDual_e
	phi_i 
	Sigma_aa_i = the variance of the ith alpha.
MEMORY:
	dp length = 2 * ka + max(
	                                        ka *n,
				ka *ka,
				2 * ka	)
*/
{
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){
		pfda_debug_msg("entering pfda_bin_cond_aa\n");
		pfda_debug_msg("w_i:\n");printmat(w_i,one,*n);
		pfda_debug_msg("ww_i:\n");printmat(ww_i,*n,*n);
		pfda_debug_msg("n: %d\n",*n);
		pfda_debug_msg("M: %d\n",*M);
		pfda_debug_msg("S1_i:\n");printmat(S1_i,*ka,one);
		pfda_debug_msg("phi_i:\n");printmat2(*n,*ka,phi_i,M);
		pfda_debug_msg("Sigma_aa_i:\n");printmat(Sigma_aa_i,*ka,*ka);
		pfda_debug_msg("ka: %d\n",*ka);
		fflush(stdout);
		}
	
	double * b = pfdaAlloc_d(*ka**ka, &dp);
	double * tmp = pfdaAlloc_d(*n**ka, &dp);
	
	dsymm_(&Left, &Upper, n, ka, &dOne, ww_i, n, phi_i, M, &dzero, tmp, n);
		if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("ww_i phi_i:\n");printmat(tmp,*n,*ka);}
	dgemm_(&Trans, &NoTrans, ka,ka,n, &dOne, phi_i, M, tmp, n, &dzero, b, ka); 
		if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("t(phi_i) ww_i phi_i:\n");printmat(b,*ka,*ka);}
	dsymm_(&Right, &Upper, ka, ka, &dOne, Sigma_aa_i, ka, b, ka, &dzero, tmp, ka);
		if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("t(phi_i) ww_i phi_i Sa_i:\n");printmat(tmp,*ka,*ka);}
	dsymm_(&Left, &Upper, ka, ka, &dOne, Sigma_aa_i, ka, tmp, ka, &dzero, aa_hat, ka);
	// pfda_matrix_inner_quadratic_form(b, phi_i, ka, M, ww_i, n, dl, dp);
	// pfda_matrix_outer_quadratic_form(aa_hat, Sigma_aa_i, ka, ka, b, ka, dl, dp);
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("aa_hat (1)\n");printmat(aa_hat,*ka,*ka);}
		
	double * a = pfdaAlloc_d(*ka,&dp);
	double * sa = pfdaAlloc_d(*ka,&dp);
	dgemv_(&Trans, n, ka, &dOne, phi_i, M, w_i, &one, &dzero, a, &one);
	dsymv_(&Upper, ka, &dOne, Sigma_aa_i, ka, a, &one, &dzero,sa, &one);
	dsyr2_(&Upper, ka, &dOne, sa, &one, S1_i, &one, aa_hat, ka);
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("aa_hat (2)\n");printmat(aa_hat,*ka,*ka);}
	
	dsyr_(&Upper, ka, &dOne, S1_i, &one, aa_hat, ka);
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("aa_hat (3)\n");printmat(aa_hat,*ka,*ka);}
	
	int ka2=*ka**ka;
	daxpy_(&ka2, &dOne, Sigma_aa_i, &one, aa_hat, &one);
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("aa_hat (4)\n");printmat(aa_hat,*ka,*ka);}
	
	pfda_fillsym(aa_hat, ka, dl);
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa)){pfda_debug_msg("aa_hat\n");printmat(aa_hat,*ka,*ka);}
	if(checkdebug(dl,debugnum_pfda_bin_cond_aa))pfda_debug_msg("exiting pfda_bin_cond_aa\n");
}

void test_pfda_bin_cond_aa(
	/* output  */
	double * const aa_hat,
	/* input */
	double const * const w_i, 
	double const * const ww_i, 
	int    const * const n,
	int    const * const M,
	double const * const S1_i,
	double const * const phi,
	double const * const Sigma_aa_i,
	int    const * const ka,	
	/* Control Values */
	int const * const phioffset,
	int const * const dl,
	double * dp){
	if(checkdebug(dl,debugnum_test_pdfa_bin_cond_aa))pfda_debug_msg("offset: %d\n",*phioffset);
	pfda_bin_cond_aa(aa_hat,w_i, ww_i, n,M,S1_i, phi+*phioffset, Sigma_aa_i, ka,dl,dp);	}

void pfda_bin_cond_bb(
	/* output  */
	double * const bb_hat,
	/* input */
	double const * const w_i, 
	double const * const ww_i, 
	int    const * const n,
	int    const * const M,
	double const * const S2_i,
	double const * const phi_i,
	double const * const Sigma_ab_i,
	double const * const Sigma_bb_i,
	int    const * const ka,
	int    const * const kb,
	
	/* Control Values */
	int const * const dl,
	double * dp)
/* DESCRIPTION
computes the  conditional expectation of \widehat{\beta\beta\trans}
INPUTS:
	w = the expected value of w given y and z.
	ww_i = the expected value of ww^t given y and z
	n = the number of observations for this subject/curve
	M = the total number of observations (number of rows of phi_i)
	S2_i = vector of length kb from the computations from pfdaDual_e
	phi_i 
	Sigma_bb_i = the variance of the ith alpha.
MEMORY:
	dp length = ka^2 +  max(
					ka * n, 
					ka * kb, 
					ka + kb	)
*/
{
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb)){
		pfda_debug_msg("entering pfda_bin_cond_bb\n");fflush(stdout);
		pfda_debug_msg("w_i:\n");printmat(w_i,one,*n);
		pfda_debug_msg("ww_i:\n");printmat(ww_i,*n,*n);
		pfda_debug_msg("n: %d\n",*n);
		pfda_debug_msg("M: %d\n",*M);
		pfda_debug_msg("S2_i:\n");printmat(S2_i,*kb,one);
		pfda_debug_msg("phi_i:\n");printmat2(*n,*ka,phi_i,M);
		pfda_debug_msg("Sigma_ab_i:\n");printmat(Sigma_ab_i,*ka,*kb);
		pfda_debug_msg("Sigma_bb_i:\n");printmat(Sigma_bb_i,*kb,*kb);
		pfda_debug_msg("ka: %d\n",*ka);
		pfda_debug_msg("kb: %d\n",*kb);
		fflush(stdout);
		}
	double * b = pfdaAlloc_d(*ka**ka, &dp);
	pfda_matrix_inner_quadratic_form(b, phi_i, ka, M, ww_i, n, dl, dp);
	pfda_matrix_inner_quadratic_form(bb_hat, Sigma_ab_i, kb, ka, b, ka, dl, dp);
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb)){pfda_debug_msg("b\n");printmat(b,*ka,*ka);
		                    pfda_debug_msg("bb_hat (1)\n");printmat(bb_hat,*ka,*ka);}
		
	double * a = pfdaAlloc_d(*ka,&dp);
	double * sa = pfdaAlloc_d(*kb,&dp);
	dgemv_(&Trans, n, ka, &dOne, phi_i, M, w_i, &one, &dzero, a, &one);
	dgemv_(&Trans, ka,kb, &dOne, Sigma_ab_i, ka, a, &one, &dzero,sa, &one);
	dsyr2_(&Upper, kb, &dOne, sa, &one, S2_i, &one, bb_hat, kb);
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb)){pfda_debug_msg("bb_hat (2)\n");printmat(bb_hat,*kb,*kb);}
	
	dsyr_(&Upper, kb, &dOne, S2_i, &one, bb_hat, kb);
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb)){pfda_debug_msg("bb_hat (3)\n");printmat(bb_hat,*kb,*kb);}
	
	int kb2=*kb**kb;
	daxpy_(&kb2, &dOne, Sigma_bb_i, &one, bb_hat, &one);
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb)){pfda_debug_msg("bb_hat (4)\n");printmat(bb_hat,*kb,*kb);}
	
	pfda_fillsym(bb_hat, kb, dl);
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb)){pfda_debug_msg("bb_hat\n");printmat(bb_hat,*kb,*kb);}
	if(checkdebug(dl,debugnum_pfda_bin_cond_bb))pfda_debug_msg("exiting pfda_bin_cond_bb\n");
}

void test_pfda_bin_cond_bb(
	/* output  */
	double * const bb_hat,
	/* input */
	double const * const w_i, 
	double const * const ww_i, 
	int    const * const n,
	int    const * const M,
	double const * const S2_i,
	double const * const phi,
	double const * const Sigma_ab_i,
	double const * const Sigma_bb_i,
	int    const * const ka,	
	int    const * const kb,	
	/* Control Values */
	int const * const phioffset,
	int const * const dl,
	double * dp){
	if(checkdebug(dl,debugnum_test_pfda_bin_cond_bb))pfda_debug_msg("offset: %d\n",*phioffset);
	pfda_bin_cond_bb(bb_hat,w_i, ww_i, n,M,S2_i, phi+*phioffset, Sigma_ab_i,Sigma_bb_i, ka, kb,dl,dp);	}

void pfda_bin_cond_ab(
	/* output  */
	double * const ab_hat,
	/* input */
	double const * const w_i, 
	double const * const ww_i, 
	int    const * const n,
	int    const * const M,
	double const * const S1_i,
	double const * const S2_i,
	double const * const phi_i,
	double const * const Sigma_aa_i,
	double const * const Sigma_ab_i,
	int    const * const ka,	
	int    const * const kb,	
	/* Control Values */
	int const * const dl,
	double * dp)
/* DESCRIPTION
computes the  conditional expectation of \widehat{\alpha\beta\trans}
INPUTS:
	w = the expected value of w given y and z.
	ww_i = the expected value of ww^t given y and z
	n = the number of observations for this subject/curve
	M = the total number of observations (number of rows of phi_i)
	S1_i = vector of length ka from the computations from pfdaDual_e
	S2_i = vector of length kb from the computations from pfdaDual_e
	phi_i 
	Sigma_ab_i = the covariance of the ith alpha and beta.
	Sigma_aa_i = the variance of the ith alpha.
MEMORY:
	dp length = k^2 + max (
	                                  ka * n,
				      ka^2 + 2 * ka + kb)
*/
{
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){
		pfda_debug_msg("entering pfda_bin_cond_ab\n");fflush(stdout);
		pfda_debug_msg("w_i:\n");printmat(w_i,one,*n);
		pfda_debug_msg("ww_i:\n");printmat(ww_i,*n,*n);
		pfda_debug_msg("n: %d\n",*n);
		pfda_debug_msg("M: %d\n",*M);
		pfda_debug_msg("S1_i:\n");printmat(S1_i,*ka,one);
		pfda_debug_msg("S2_i:\n");printmat(S2_i,*kb,one);
		pfda_debug_msg("phi_i:\n");printmat2(*n,*ka,phi_i,M);
		pfda_debug_msg("Sigma_aa_i:\n");printmat(Sigma_aa_i,*ka,*ka);
		pfda_debug_msg("Sigma_ab_i:\n");printmat(Sigma_ab_i,*ka,*kb);
		pfda_debug_msg("ka: %d\n",*ka);
		pfda_debug_msg("kb: %d\n",*kb);
		fflush(stdout);
		}

	double * b = pfdaAlloc_d(*ka**ka, &dp);
	pfda_matrix_inner_quadratic_form(b, phi_i, ka, M, ww_i, n, dl, dp);
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("b\n");printmat(b,*ka,*ka);}
	double * c = pfdaAlloc_d(*ka**ka,&dp);
	dsymm_(&Left, &Upper, ka, ka, &dOne, Sigma_aa_i, ka, b, ka, &dzero,c, ka);
	dgemm_(&NoTrans, &NoTrans, ka, kb, ka, &dOne, c, ka, Sigma_ab_i, ka, &dOne, ab_hat, ka); 
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("ab_hat (1)\n");printmat(ab_hat,*ka,*kb);}
	
	double * a =pfdaAlloc_d(*ka, &dp);
	dgemv_(&Trans, n, ka, &dOne, phi_i, M, w_i, &one, &dzero, a, &one);
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("ab_hat (2)\n");printmat(ab_hat,*ka,*kb);}
	
	double * sb = pfdaAlloc_d(*kb,&dp);
	dgemv_(&Trans, ka, kb, &dOne, Sigma_ab_i, ka, a, &one, &dzero, sb, &one);
	dger_(ka, kb, &dOne, S1_i, &one, sb, &one, ab_hat, ka);
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("ab_hat (3)\n");printmat(ab_hat,*ka,*kb);}

	double * sa = pfdaAlloc_d(*ka,&dp);
	dsymv_(&Upper, ka, &dOne, Sigma_aa_i, ka, a, &one, &dzero,sa, &one);
	dger_(ka, kb, &dOne, sa, &one, S2_i, &one, ab_hat, ka);
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("ab_hat (4)\n");printmat(ab_hat,*ka,*kb);}
	
	dger_(ka, kb, &dOne, S1_i, &one, S2_i, &one, ab_hat, ka);
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("ab_hat (5)\n");printmat(ab_hat,*ka,*kb);}
	
	int kab=*ka**kb;
	daxpy_(&kab, &dOne, Sigma_ab_i, &one, ab_hat, &one);
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab)){pfda_debug_msg("ab_hat (final)\n");printmat(ab_hat,*ka,*kb);
		pfda_debug_msg("leaving pfda_bin_cond_ab");}
}    
void test_pfda_bin_cond_ab(
	/* output  */
	double * const ab_hat,
	/* input */
	double const * const w_i, 
	double const * const ww_i, 
	int    const * const n,
	int    const * const M,
	double const * const S1_i,
	double const * const S2_i,
	double const * const phi,
	double const * const Sigma_aa_i,
	double const * const Sigma_ab_i,
	int    const * const ka,	
	int    const * const kb,	
	/* Control Values */
	int const * const phioffset,
	int const * const dl,
	double * dp){
	if(checkdebug(dl,debugnum_pfda_bin_cond_ab))pfda_debug_msg("offset: %d\n",*phioffset);
	pfda_bin_cond_ab(ab_hat, w_i, ww_i, n, M, S1_i, S2_i, phi+*phioffset, Sigma_aa_i, Sigma_ab_i, ka, kb, dl, dp);	}


