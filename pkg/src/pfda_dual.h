/*! 	\file pfda_dual.h
	\brief Algorithm Steps for the Bivariate Linear Algorithm
*/

/*! 	\defgroup dual Dual Algorithm Steps.
	steps for pfda_dual()
*/

#ifndef PFDA_DUAL_H_
#define PFDA_DUAL_H_
extern void pfdaDual_e(
	/* input Values */
	double const * const y, 
	double const * const z, 
	int    const * const nobs,
	int    const * const M, 
	int    const * const N, 
	int    const * const ka, 
	int    const * const kb, 
	double const * const B, 
	int    const * const p,
	/* State Values */
	double const * const tm, 
	double const * const tn, 
	double const * const tf, 
	double const * const tg, 
	double const * const Da, 
	double const * const Db, 
	double const * const Lambda,
	double const * const seps, 
	double const * const sxi, 
	double * const Alpha, 
	double * const Beta, 
	double * const Sigma_aa, 
	double * const Sigma_ab, 
	double * const Sigma_bb,
	/* Control Values */
	int const * const dl, double * dp, int * ip);
extern void pfda_cond_ab(
	double * const ab_hat,
	int const * const ld,
	double const * const Alpha,
	int const * const ka,
	double const * const Beta, 
	int const * const kb,
	double const * const Sigma_ab,
	int const * const dl
	);
extern void pfda_sum_cond_ab(
	double * const sum_ab_hat,
	int const * const N,
	double const * const Alpha,
	int const * const ka,
	double const * const Beta, 
	int const * const kb,
	double const * const Sigma_ab,
	int const * const dl
	);	
void dual_gen_sigmas(
	      double * const Saa,
	      double * const Sab,
	      double * const Sbb,
	const double * const phi,
	const double * const psi,
	const double * const lambda,
	const double * const Da,
	const double * const Db,
	const double * const sep,
	const double * const sxi,
	const int    * const M,
	const int    * const ni,
	const int    * const ka,
	const int    * const kb,
	const int * const dl,
	double*dp, int * ip);
void pfdaDual_m5_2(
	double       * const Lambda,
	double       * const trans_f, //will hold the chol decomp after run.
	int    const * const ka,
	double const * const trans_g,
	int    const * const kb,
	int const * const dl, double*dp,int*ip);
void dual_bc_5(
	double * lambda,
	const double * const aa,
	const double * const ab,
	const int * const N,
	const int * const ka,
	const int * const kb,
	const int * const dl,double * dp, int * ip);
void pfdaDual_m5_2(
	double       * const Lambda,
	double       * const trans_f, //will hold the chol decomp after run.
	int    const * const ka,
	double const * const trans_g,
	int    const * const kb,
	int const * const dl, double*dp,int*ip);
#endif /* PFDA_DUAL_H_ */
