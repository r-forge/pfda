/*!	\file pfda_utils.h
	\brief general utility functions for pfda package
*/
#ifndef PFDA_UTILS_H_
#define PFDA_UTILS_H_

#include "pfda_const.h"

#define MMIN(x,y) (((x)<(y))?(x):(y))
#define MMAX(x,y) (((x)>(y))?(x):(y))

/*! 	\defgroup utils Utility Functions */
/*\@{*/
extern double* Idiag(int i);
extern void resetI(double* A,const int n);
extern int local_dspev(const int* n, double* A, double* W, double* Z, const int *ldz, double* work);
extern void positivefirstrow(double *A,const int nrow,const int ncol);
extern void pfda_transpose(double * const A, const int nrow, const int ncol,int const * const dl, double * const dp);
extern void myinner(double* A,int* ncol,int* nrow,double* C,int* ldc);
extern double tr(double * A,int n);
extern void setMinVar(double* const D, int const*const k, double const*const var);
extern void _dqfmm(const char* UPLO, const int* m,const int* k,const double* alpha, const double* A, const int* lda, const double* S, const int* lds, const double* beta, double *Q, const int* ldq, const int debug);
extern void _dtqimm(const char* UPLO, const int* m,const int* k,const double* alpha, const double* A, const int* lda, double* S, const int* lds, const double* beta, double *Q, const int* ldq, int *info);
extern void _dqimm(const char* UPLO, const int* m,const int* k,const double* alpha, const double* A, const int* lda, double* S, const int* lds, const double* beta, double *Q, const int* ldq, int *info);
extern void reverse(const int *N, double *V);
extern void reverseColumns(const int *ncol, double *A, const int* lda);
extern void addDiag(int const * const N, double const * const D, double * const A, int const *const lda);
extern void addInvDiag(int const * const N, double const * const D, double * const A, int const *const lda);
extern void _dchsm(const char *side, const char *uplo, const int *M, const int *N,const double *alpha, const double *A,const int *lda, double *B, const int *ldb);
extern void pfda_fillsym(double * const sym,int const * const n, int const * const dl);
extern void pfda_sym_inverse(double * const A,int const * const n, int * sr, int const * const dl, double * dp, int * ip);
int pfda_arrMAXi(int * arr, int * length);
extern void pfda_matrix_outer_quadratic_form(
	double * const Q,
	double const * const X,
	int const * const nrowx,
	int const * const ldx,
	double const * const S,
	int const * const nrows,
	int const * const dl,
	double * dp
	);
extern void pfda_matrix_inner_quadratic_form(
	double       * const Q,
	double const * const X,	int    const * const ncolx,	int    const * const ldx,
	double const * const S,	int    const * const nrows,
	int const * const dl,
	double * dp);
extern void pfda_eigens(
	double * const S,	int const * const ns,
	double * const vectors,	double * const values,	int const * const k,
	int const * const dl, double * dp, int * ip);
extern void pfda_uppersymadd(
	int const * const n,
	double * const x, int const*const ldx,
	double const * const alpha,
	double const * const y, int const*const ldy,
	int const * const dl);
extern void pfda_computebtb(
	double *btb,
	int const * const N,
	double const * const B,
	int const * const M,
	int const * const p,
	int const * const nobs,
	int const * const dl);
extern void pfda_computeResid(
	double       * const Ry,
	double const * const y, 		
	int    const * const nobs,
	int    const * const M, 
	int    const * const N, 
	int    const * const k, 
	double const * const B, 
	int    const * const p,
	double const * const tm, 
	double const * const tf, 
	double const * const Alpha, 
	int const * const dl, double*dp);
extern void pfda_rtruncnorm(
	double       * const x,
	int    const * const n,
	double const * const lower, int    const * const nl,
	double const * const upper,	int    const * const nu,
	double const * const mean ,	int    const * const nm,
	double const * const sd   ,	int    const * const ns,
	int const * const dl);
void pfda_gen_truncnorm(
	double * x,
	int const * const n,
	double const * const c,
	int const * const dl);
void gen_orthog(
	double * const tf,
	double * const alpha,
	double * const Da,
	double * const trans,
	const double * const avgaa,
	const int * const N,
	const int * const ka,
	const int * const p,
	const int * const dl,
	double * dp, int * ip);

/*\@}*/
#endif /* PFDA_UTILS_H_ */
