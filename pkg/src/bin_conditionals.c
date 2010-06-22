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
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa){
		pfda_debug_msg("entering pfda_bin_cond_aa\n");
		pfda_debug_argvec(w_i,n);
		pfda_debug_argmat(ww_i,*n,*n);
		pfda_debug_arg(*n);
		pfda_debug_arg(*M);
		pfda_debug_argmat(S1_i,*ka,one);
		pfda_debug_argmat2(*n,*ka,phi_i,M);
		pfda_debug_argmat(Sigma_aa_i,*ka,*ka);
		pfda_debug_arg(*ka);
		fflush(stdout);
		}
	
	double * b = pfdaAlloc_d(*ka**ka, &dp);
	double * tmp = pfdaAlloc_d(*n**ka, &dp);
	
	dsymm_(&Left, &Upper, n, ka, &dOne, ww_i, n, phi_i, M, &dzero, tmp, n);
		pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(tmp,*n,*ka);}
	dgemm_(&Trans, &NoTrans, ka,ka,n, &dOne, phi_i, M, tmp, n, &dzero, b, ka); 
		pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(b,*ka,*ka);}
	dsymm_(&Right, &Upper, ka, ka, &dOne, Sigma_aa_i, ka, b, ka, &dzero, tmp, ka);
		pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(tmp,*ka,*ka);}
	dsymm_(&Left, &Upper, ka, ka, &dOne, Sigma_aa_i, ka, tmp, ka, &dzero, aa_hat, ka);
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(aa_hat,*ka,*ka);}
	
	double * a = pfdaAlloc_d(*ka,&dp);
	double * sa = pfdaAlloc_d(*ka,&dp);
	dgemv_(&Trans, n, ka, &dOne, phi_i, M, w_i, &one, &dzero, a, &one);
	dsymv_(&Upper, ka, &dOne, Sigma_aa_i, ka, a, &one, &dzero,sa, &one);
	dsyr2_(&Upper, ka, &dOne, sa, &one, S1_i, &one, aa_hat, ka);
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(aa_hat,*ka,*ka);}
	
	dsyr_(&Upper, ka, &dOne, S1_i, &one, aa_hat, ka);
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(aa_hat,*ka,*ka);}
	
	int ka2=*ka**ka;
	daxpy_(&ka2, &dOne, Sigma_aa_i, &one, aa_hat, &one);
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(aa_hat,*ka,*ka);}
	
	pfda_fillsym(aa_hat, ka, dl);
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa){pfda_debug_argmat(aa_hat,*ka,*ka);}
	pfda_debug_cdl(debugnum_pfda_bin_cond_aa)pfda_debug_msg("exiting pfda_bin_cond_aa\n");
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
	pfda_debug_cdl(debugnum_test_pdfa_bin_cond_aa)pfda_debug_arg(*phioffset);
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
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb){
		pfda_debug_msg("entering pfda_bin_cond_bb\n");fflush(stdout);
		pfda_debug_argmat(w_i,one,*n);
		pfda_debug_argmat(ww_i,*n,*n);
		pfda_debug_arg(*n);
		pfda_debug_arg(*M);
		pfda_debug_argmat(S2_i,*kb,one);
		pfda_debug_argmat2(*n,*ka,phi_i,M);
		pfda_debug_argmat(Sigma_ab_i,*ka,*kb);
		pfda_debug_argmat(Sigma_bb_i,*kb,*kb);
		pfda_debug_arg(*ka);
		pfda_debug_arg(*kb);
		fflush(stdout);
		}
	double * b = pfdaAlloc_d(*ka**ka, &dp);
	pfda_matrix_inner_quadratic_form(b, phi_i, ka, M, ww_i, n, dl, dp);
	pfda_matrix_inner_quadratic_form(bb_hat, Sigma_ab_i, kb, ka, b, ka, dl, dp);
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb){pfda_debug_argmat(b,*ka,*ka);
		                    pfda_debug_argmat(bb_hat,*ka,*ka);}
		
	double * a = pfdaAlloc_d(*ka,&dp);
	double * sa = pfdaAlloc_d(*kb,&dp);
	dgemv_(&Trans, n, ka, &dOne, phi_i, M, w_i, &one, &dzero, a, &one);
	dgemv_(&Trans, ka,kb, &dOne, Sigma_ab_i, ka, a, &one, &dzero,sa, &one);
	dsyr2_(&Upper, kb, &dOne, sa, &one, S2_i, &one, bb_hat, kb);
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb){pfda_debug_argmat(bb_hat,*kb,*kb);}
	
	dsyr_(&Upper, kb, &dOne, S2_i, &one, bb_hat, kb);
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb){pfda_debug_argmat(bb_hat,*kb,*kb);}
	
	int kb2=*kb**kb;
	daxpy_(&kb2, &dOne, Sigma_bb_i, &one, bb_hat, &one);
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb){pfda_debug_argmat(bb_hat,*kb,*kb);}
	
	pfda_fillsym(bb_hat, kb, dl);
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb){pfda_debug_argmat(bb_hat,*kb,*kb);}
	pfda_debug_cdl(debugnum_pfda_bin_cond_bb)pfda_debug_msg("exiting pfda_bin_cond_bb\n");
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
	pfda_debug_cdl(debugnum_test_pfda_bin_cond_bb)pfda_debug_arg(*phioffset);
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
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){
		pfda_debug_msg("entering pfda_bin_cond_ab\n");fflush(stdout);
		pfda_debug_argmat(w_i,one,*n);
		pfda_debug_argmat(ww_i,*n,*n);
		pfda_debug_arg(*n);
		pfda_debug_arg(*M);
		pfda_debug_argmat(S1_i,*ka,one);
		pfda_debug_argmat(S2_i,*kb,one);
		pfda_debug_argmat2(*n,*ka,phi_i,M);
		pfda_debug_argmat(Sigma_aa_i,*ka,*ka);
		pfda_debug_argmat(Sigma_ab_i,*ka,*kb);
		pfda_debug_arg(*ka);
		pfda_debug_arg(*kb);
		fflush(stdout);
		}

	double * b = pfdaAlloc_d(*ka**ka, &dp);
	pfda_matrix_inner_quadratic_form(b, phi_i, ka, M, ww_i, n, dl, dp);
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(b,*ka,*ka);}
	double * c = pfdaAlloc_d(*ka**ka,&dp);
	dsymm_(&Left, &Upper, ka, ka, &dOne, Sigma_aa_i, ka, b, ka, &dzero,c, ka);
	dgemm_(&NoTrans, &NoTrans, ka, kb, ka, &dOne, c, ka, Sigma_ab_i, ka, &dOne, ab_hat, ka); 
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(ab_hat,*ka,*kb);}
	
	double * a =pfdaAlloc_d(*ka, &dp);
	dgemv_(&Trans, n, ka, &dOne, phi_i, M, w_i, &one, &dzero, a, &one);
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(ab_hat,*ka,*kb);}
	
	double * sb = pfdaAlloc_d(*kb,&dp);
	dgemv_(&Trans, ka, kb, &dOne, Sigma_ab_i, ka, a, &one, &dzero, sb, &one);
	dger_(ka, kb, &dOne, S1_i, &one, sb, &one, ab_hat, ka);
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(ab_hat,*ka,*kb);}

	double * sa = pfdaAlloc_d(*ka,&dp);
	dsymv_(&Upper, ka, &dOne, Sigma_aa_i, ka, a, &one, &dzero,sa, &one);
	dger_(ka, kb, &dOne, sa, &one, S2_i, &one, ab_hat, ka);
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(ab_hat,*ka,*kb);}
	
	dger_(ka, kb, &dOne, S1_i, &one, S2_i, &one, ab_hat, ka);
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(ab_hat,*ka,*kb);}
	
	int kab=*ka**kb;
	daxpy_(&kab, &dOne, Sigma_ab_i, &one, ab_hat, &one);
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab){pfda_debug_argmat(ab_hat,*ka,*kb);
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
	pfda_debug_cdl(debugnum_pfda_bin_cond_ab)pfda_debug_arg(*phioffset);
	pfda_bin_cond_ab(ab_hat, w_i, ww_i, n, M, S1_i, S2_i, phi+*phioffset, Sigma_aa_i, Sigma_ab_i, ka, kb, dl, dp);	}


