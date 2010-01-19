/*! 	\file pfda_single.h
	\brief Algorithm Steps for the Linear Univeriate Algorithm
*/

/*! 	\defgroup single Single Algorithm Steps
	steps for pfdaSingle()
	also holds some conditionals computations.
*/

#ifndef PFDA_SINGLE_H_
#define PFDA_SINGLE_H_
extern void pfda_s_i(
	double       * const tm,
	double       * const tf,
	double       * const alpha,
	double       * const Da,
	double       * const aa,
	double       * const sigma,
	double const * const y, 
	int    const * const nobs,
	int    const * const M, 
	int    const * const N, 
	int    const * const k, 
	double const * const B, 
	double const * const btb,
	int    const * const p, 
	double const * const minV,
	int const * const dl,
	double * dp, int * ip);
extern void pfdaSingle_e_1(
	double       * const Alpha_i,
	int    const * const ldalpha,
	double       * const Sigma_i,
	double const * const Ryi,
	int    const * const nobs_i,
	int    const * const M,
	int    const * const ka,
	double const * const phii,
	double const * const Da,
	double const * const inv_epsilon,
	int const * const dl, double *dp , int * ip);
extern void pfdaSingle_e(
	double       * const Alpha,
	double       * const Sigma_aa,
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
	double const * const seps,
	const int * const dl, double *dp , int * ip);
extern void pfdaSingle_m5(
	int    const * const N,
	double       * const delta,
	int    const * const k,
	double const * const Sigma_dd,
	double       * const Theta,
	double       * const D,
	int    const * const p,
	double const * const minimum_variance,
	int const * const dl, double * dp,int * ip
	);

#endif /* PFDA_SINGLE_H_ */
