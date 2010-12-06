/*! 	\file pfda_gen.h
	\brief Algorithm Steps use Generally
*/

/*! 	\defgroup generalsteps General Steps.
	These steps are used by at least two of the general algorithms (may be nested deeper, and have wrapper functions)
*/

#ifndef PFDA_GEN_H_
#define PFDA_GEN_H_
extern void pfdafindinits(double const * const y, double const * const z, int const * const nobs,int const * const M, int const * const N, int const * const ka,int const * const kb,double const * const B, int const * const p, double const * const delta, 
double* const tm,double* const tn,double* const tf, double* const tg, double* const Da, double* const Db, double* const L, double* const seps, double* const sxi, double* const Alpha, double* const Beta, double* const Saa, double* const Sab, double* const Sbb); 

extern void pfdaFindInitSingle(
	double const * const t, 
	double const * const y, 
	int const * const nobs,
	int const * const M, 
	int const * const N, 
	int const * const ka, 
	double const * const B, 
	int const * const p, 
	double const * const delta, 
	double * const tm, 
	double * const tf, 
	double * const Da, 
	double * const seps, 
	double * const Alpha,
	const int *const dl); 

extern void pfda_m1(
	double       * const seps, 
	/* Inputs Values */
	double const * const y, 
	int    const * const nobs,
	int    const * const M, 
	int    const * const N, 
	int    const * const ka, 
	double const * const B, 
	int    const * const p,
	double const * const delta, 
	double const * const tm, 
	double const * const tf, 
	double const * const Alpha, 
	double const * const Sigma_aa, 
	/* Auxilliary */
	int const * const dl, double * dp
	);
extern void pfda_m2(
	double * const tm, 
	/* input Values */
	double const * const y, 
	const int * const  nobs,
	const int * const  M, 
	const int * const  N, 
	const int * const  ka, 
	double const * const B, 
	const int * const  p,
	double const * const lm, 
	double const * const ISD,
	double const * const tf, 
	double const * const seps, 
	double const * const Alpha, 
	int const * const dl, double * dp);
extern void pfda_m3_core(
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
	double const * const ISD,
	double const * const tm,
	double const * const seps,
	double const * const Alpha,
	double const * const btb,
	int const * const dl, double * dp, int * ip);
extern void pfda_m3(
	double * const tf, 
	/* input Values */
	double const * const y, 
	const int * const  nobs,
	const int * const  M, 
	const int * const  N, 
	const int * const  ka, 
	double const * const B, 
	int    const * const p,
	double const * const lf, 
	double const * const ISD,
	double const * const tm, 
	double const * const seps, 
	double const * const Alpha, 
	double const * const Sigma_aa,
	double const * const btb,
	int const * const dl, double * dp, int * ip);
extern void pfda_m5_0(
	double       * const delta,
	double const * const sum_dd_hat,
	int    const * const N,
	int    const * const k,
	double const * const Sigma_dd,
	double       * const Theta,
	double       * const D,
	int    const * const p,
	double       * const Transformation,
	double const * const minimum_variance,
	int const * const dl, double * dp, int * ip);
extern void pfda_m5_1(
	int const * const N,
	double * const delta,
	int const * const k,
	double const * const Sigma_dd,
	double * const Theta,
	double * const D,
	int const * const p,
	double * const Transformation,
	double const * const minimum_variance,
	int const * const dl, double * dp, int * ip	);
extern void pfda_cond_dd(
	double * const dd_hat, 
	double const * const Delta_i, 
	int const * const N,
	double const * const Sigma_dd_i,
	int const * const kd,
	int const * const dl
	);
extern void pfda_sum_cond_dd(
	double * const sum_dd_hat, 
	int const * const N,
	double const * const Delta, 
	double const * const Sigma_dd,
	int const * const kd,
	int const * const dl
	);
extern void pfda_bin_cond_aa(
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
	double * dp);
#endif /* PFDA_GEN_H_ */
