/*! \file conditionals.c
//  by Andrew Redd
// Part of paired functional data analysis through principle component analysis.
// algorithms for finding the conditional expectations of alpha and beta.
*/

#include "pfda.h"

void pfda_cond_dd(
	double * const dd_hat, 
	double const * const Delta_i, 
	int const * const lddelta,
	double const * const Sigma_dd_i,
	int const * const kd,
	int const * const dl
	)
/* DESCRIPTION
Generic function for either Alpha or Beta conditional given Y and/or Z,  Returns a full matrix not only a symatric one but could be use with symatric ones. A little less efficient but more easily tractable.
INPUTS:
	dd_hat (ka×ka) - output variable, assumed to hold zeros, but not enforced to hold zeros, so could be use to hold a sum tally.
	Delta _i (kd) - represents either Alpha or Beta matrix of principle component scores for Y/Z.  ROW i represents delta_i
	N - the number of subjects. Also the number of rows of Delta.
	Sigma_dd _i (kd×kd) - Three dimentional array.  each ith block of length ka^2 represents the variance for alpha_i
	kd  - number of columns of Delta, ie. the number of principle components for Y/Z.
	dl - integer array for printing debugging information.
OUTPUT:
	dd_hat
	\[ %tex description
		\widehat{\alpha_i\alpha_i\trans} = E(\alpha_i\alpha_i\trans|Y_i,Z_i) = \hat\alpha_i\hat\alpha_i\trans + \Simga_{i,\alpha\alpha}
	\]
debug number:	301
*/
{
	pfda_debug_cdl(debugnum_cond_dd){
		pfda_debug_msg("pfda_cond_dd:\n");
		pfda_debug_msg("dd_hat:\n"); printmat(dd_hat, *kd, *kd);
		pfda_debug_msg("Delta:\n"); printmat2(one, *kd, Delta_i, lddelta);
		pfda_debug_msg("lddelta:\n%d\n\n",*lddelta );
		pfda_debug_msg("Sigma_dd[i]:\n");  printmat(Sigma_dd_i, *kd, *kd);
		pfda_debug_msg("kd:\n%d\n\n",*kd);
	}
	dsyr_(&Upper, kd, &dOne, Delta_i, lddelta, dd_hat,  kd);
	pfda_fillsym(dd_hat, kd, dl);
	int const kd2 = *kd**kd;
	daxpy_(&kd2, &dOne, Sigma_dd_i, &one, dd_hat, &one); // dd_hat += Sigma_dd[i]
}
void test_pfda_cond_dd(	
	double * const dd_hat, 
	int const * const i, 
	int const * const N,
	double const * const Delta, 
	double const * const Sigma_dd,
	int const * const kd,
	int const * const dl
) {pfda_cond_dd(dd_hat, Delta+*i, N, Sigma_dd+*i**kd**kd, kd,dl);}
void pfda_cond_ab(
	double * const ab_hat,
	int const * const ld,
	double const * const Alpha_i,
	int const * const ka,
	double const * const Beta_i, 
	int const * const kb,
	double const * const Sigma_ab_i,
	int const * const dl
	)
/* DESCRIPTION
Generic function for either Alpha or Beta conditional given Y and/or Z,  Returns a full matrix not only a symatric one but could be use with symatric ones. A little less efficient but more easily tractable.
INPUTS:
	ab_hat (ka×ka) - output variable, assumed to hold zeros, but not enforced to hold zeros, so could be use to hold a sum tally.
	N - the number of subjects. Also the number of rows of Delta.
	i - the index to retrieve
	Alpha (N×ka) - matrix of principle component scores for Y.  ROW i represents alpha_i
	ka  - number of columns of Alpha, ie. the number of principle components for Y.
	Beta (N×ka) - matrix of principle component scores for Z.  ROW i represents beta_i
	kb  - number of columns of Beta, ie. the number of principle components for Z.
	Sigma_ab (N×ka×kb) - Three dimentional array.  each ith block of length ka*kb represents the covariance for alpha_i and beta_i
	dl - integer array for printing debugging information.
OUTPUT:
	dd_hat
	\[ %tex description
		\widehat{\alpha_i\beta_i\trans} = E(\alpha_i\beta_i\trans|Y_i,Z_i) = \hat\alpha_i\hat\beta_i\trans + \Simga_{i,\alpha\beta}
	\]
dl:	302
*/
{
	pfda_debug_cdl(debugnum_cond_ab){
		pfda_debug_msg("pfda_cond_ab:\n");
		pfda_debug_msg("ab_hat:\n"); printmat(ab_hat, *ka, *kb);
		pfda_debug_msg("N:\n%d\n\n",*ld );
		pfda_debug_msg("Alpha:\n"); printmat2(one, *ka, Alpha_i, ld);
		pfda_debug_msg("ka:\n%d\n\n",*ka);
		pfda_debug_msg("Beta:\n");  printmat2(one, *kb, Beta_i , ld);
		pfda_debug_msg("kb:\n%d\n\n",*kb);
		pfda_debug_msg("Sigma_ab[i]:\n");  printmat(Sigma_ab_i, *ka, *kb);
	}
	int const kab = *ka**kb;
	dger_(ka, kb, &dOne, Alpha_i, ld, Beta_i, ld, ab_hat, ka);  //dd_hat += delta_i delta_i\trans
	daxpy_(&kab, &dOne, Sigma_ab_i, &one, ab_hat, &one); // dd_hat += Sigma_dd[i]
}
void test_pfda_cond_ab(
	double * const ab_hat,
	int const * const i,
	int const * const N,
	double const * const Alpha,
	int const * const ka,
	double const * const Beta, 
	int const * const kb,
	double const * const Sigma_ab,
	int const * const dl
	) {pfda_cond_ab(ab_hat, N, Alpha+*i, ka, Beta+*i, kb, Sigma_ab+*i**ka**kb, dl);}

void pfda_sum_cond_dd(
	double * const sum_dd_hat, 
	int const * const N,
	double const * const Delta, 
	double const * const Sigma_dd,
	int const * const kd,
	int const * const dl
	)
{
	/* DESCRIPTION
	Simple wrapper function for convienience.
	Generic function for either Alpha or Beta conditional given Y and/or Z,  Returns a full matrix not only a symatric one but could be use with symatric ones. A little less efficient but more easily tractable.	INPUTS:
		dd_hat (ka×ka) - output variable, assumed to hold zeros, but not enforced to hold zeros, so could be use to hold a sum tally.
		N - the number of subjects. Also the number of rows of Delta.
		Delta (N×kd) - represents either Alpha or Beta matrix of principle component scores for Y/Z.  ROW i represents delta_i
		Sigma_dd (N×kd×kd) - Three dimentional array.  each ith block of length ka^2 represents the variance for alpha_i
		kd  - number of columns of Delta, ie. the number of principle components for Y/Z.
		dl - integer array for printing debugging information.
	OUTPUT:
		dd_hat
		\[ %tex description
			\sum\limits_{i=1}^n\widehat{\alpha_i\alpha_i\trans} = \sum\limits_{i=1}^nE(\alpha_i\alpha_i\trans|Y_i,Z_i) = \sum\limits_{i=1}^n\hat\alpha_i\hat\alpha_i\trans + \Simga_{i,\alpha\alpha}
		\]
	dl:	303
	*/
	int const kd2 = *kd**kd;
	for(int obsnumber=0;obsnumber<*N;obsnumber++){
		pfda_cond_dd(sum_dd_hat, Delta+obsnumber, N, Sigma_dd+obsnumber*kd2, kd, dl);
	}
}
void pfda_sum_cond_ab(
	double * const sum_ab_hat,
	int const * const N,
	double const * const Alpha,
	int const * const ka,
	double const * const Beta, 
	int const * const kb,
	double const * const Sigma_ab,
	int const * const dl
	)
/* DESCRIPTION
simple wrapper function for the sum.
Generic function for either Alpha or Beta conditional given Y and/or Z,  Returns a full matrix not only a symatric one but could be use with symatric ones. A little less efficient but more easily tractable.
INPUTS:
	ab_hat (ka×ka) - output variable, assumed to hold zeros, but not enforced to hold zeros, so could be use to hold a sum tally.
	N - the number of subjects. Also the number of rows of Delta.
	Alpha (N×ka) - matrix of principle component scores for Y.  ROW i represents alpha_i
	ka  - number of columns of Alpha, ie. the number of principle components for Y.
	Beta (N×ka) - matrix of principle component scores for Z.  ROW i represents beta_i
	kb  - number of columns of Beta, ie. the number of principle components for Z.
	Sigma_ab (N×ka×kb) - Three dimentional array.  each ith block of length ka*kb represents the covariance for alpha_i and beta_i
	dl - integer array for printing debugging information.
OUTPUT:
	dd_hat
	\[ %tex description
		\sum\limits_{i=1}^n\widehat{\alpha_i\beta_i\trans} = \sum\limits_{i=1}^nE(\alpha_i\beta_i\trans|Y_i,Z_i) = \sum\limits_{i=1}^n\hat\alpha_i\hat\beta_i\trans + \Simga_{i,\alpha\beta}
	\]
dl:	302
*/
{	
	int obsnumber;
	for(obsnumber=0;obsnumber<*N;obsnumber++){
		pfda_cond_ab(sum_ab_hat, N, Alpha+obsnumber, ka, Beta+obsnumber, kb, Sigma_ab+obsnumber**ka**kb, dl);
	}
}
