/*! 	\file pfda_bin_single.h
	\brief Algorithm Steps for the Univariate Binary Algorithm
*/

#ifndef PFDA_BIN_SINGLE_H_
#define PFDA_BIN_SINGLE_H_
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
	double * dp, int * ip);
#endif
