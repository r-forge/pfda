//calcium.c
#include "pfda_error.h"
#include "pfda_gen.h"
#include "pfda_memory.h"
#include "pfda_single.h"
#include "pfda_utils.h"
#include "pfda_dual.h"
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

void zero_lambda(
	double       * const lambda,
	int const * const kg,
	int const * const kd){
	for(int i=0;i<*kg**kd;i++)lambda[i]=0.0;
}

/*! Computes residuals for the calcium model
 *
 *	@MEMORY
 *		- dp = M*max(kg,kd)
 */
void dual_ca_resid(
	double       * const Ry,
	double const * const y,
	double const * const Z,
	double const * const Bt,
	double const * const Bx,
	double const * const tz,
	double const * const tt,
	double const * const tx,
	double const * const tf,
	double const * const tg,
	double const * const gamma,
	double const * const delta,
	int    const * const nobs,
	int    const * const N,
	int    const * const M,
	int    const * const kz,
	int    const * const kd,
	int    const * const kg,
	int    const * const pt,
	int    const * const px,
	int const * const dl, double*dp)
{
	if(Ry!=y)dcopy_(M,y,&one,Ry,&one);
	if(Z && *kz)dgemv_(&NoTrans, M, kz, &mOne, Z, M, tz, &one, &dOne, Ry, &one);
	pfda_computeResid(Ry, Ry, nobs, M, N, kg, Bt, pt, tt, tf, gamma, dl, dp);
	pfda_computeResid(Ry, Ry, nobs, M, N, kd, Bx, px, tx, tg, delta, dl, dp);
}

/*! Makes vector v orthogonal to all columns in matrix U

*/
void u_orthogonalize(
	double * const v,        /// \param v the vector to be orthogonalized
	const double * const U,  /// \param U is assumed orthogonal, NULL pointer is allowed.
	const int * const nrow,  /// \param nrow is the number of rows of U and the length of v.
	const int * const ncol   /// \param ncol refers to the number of columns of U, NULL pointers and zero values are allowed
	)
{
	if(U && ncol && *ncol)for(int i=0; i<*ncol;i++){
		double a = -1.0*ddot_(nrow, v, &one, U+i**nrow, &one);
		daxpy_(nrow, &a, U+i**nrow, &one, v, &one);		
	}
	double b = 1.0/dnrm2_(nrow, v, &one);
	dscal_(nrow,&b, v, &one);	
}


/*!  Solves the matrix in a block format.  A, and C are symetric, B is the off-diagonal
@MEMORY
	- dp  = a^2 + c^2 + max(a,c)*10 + 2*a*c
	- ip = max(a,c)
*/
void gen_symblock_solve(
	double * A,
	double * B,
	double * C,
	const int * const a,
	const int * const c,
	double * dp, int * ip
)
{
	int a2 = *a**a;
	int sr = 0, lwork = *a*10;
	double * Ai = pfdaAlloc_d(a2,&dp);
	resetI(Ai,*a);
	dsysv_(&Upper, a, a, A, a, ip, Ai, a, dp, &lwork, &sr);
	// dcopy_(&a2, A, &one, Ai, &one);
	// pfda_sym_inverse(Ai, a, &sr, NULL, dp, ip);

	int cc = *c**c;
	double * BAiB = pfdaAlloc_d(*c**c, &dp);
	double * AiB = pfdaAlloc_d(*a**c,&dp);
	pfda_matrix_inner_quadratic_form(BAiB, B, c, a, Ai, a, NULL, AiB);
	dscal_(&cc, &mOne, BAiB, &one);
	daxpy_(&cc, &dOne, C, &one, BAiB, &one);
	lwork=10**c;
	resetI(C,*c);
	dsysv_(&Upper, c, c, BAiB, c, ip, C, c, dp, &lwork, &sr);

	dsymm_(&Right,&Upper, a,c, &mOne, C,c, AiB, a, &dzero, B, a);

	pfda_matrix_outer_quadratic_form(A, AiB, a, a, C, c, NULL, dp);
	daxpy_(&a2, &dOne, Ai, &one, A, &one);
	pfda_fillsym(A,a, NULL);
}

/*! computes the initial values for the dual principle component single response model

@MEMORY
	- dp = pt*pt*N + px*px*N + 2*p^2 + M + M*k+ p*N + 8*p
	- ip = kz
*/
void dual_ca_i(
	      double * const tz,
	      double * const tt,
	      double * const tx,
	      double * const tf,
	      double * const tg,
	      double * const gamma,
	      double * const delta,
	      double * const lambda,
	      double * const Dg,
	      double * const Dd,
	      double * const sigma,
	      double * const gg,
	      double * const gd,
	      double * const dd,
	const double * const y,
	const double * const z,
	const double * const Bt,
	const double * const Bx,
	const int    * const nobs,
	const int    * const N,
	const int    * const M,
	const int    * const kz,
	const int    * const kg,
	const int    * const kd,
	const int    * const pt,
	const int    * const px,
	const double * const minV,
	const int * const dl,
	double * dp, int * ip)
{
	double * btbt = pfdaAlloc_d(*pt**pt**N, &dp);
	pfda_computebtb(btbt, N, Bt, M, pt, nobs, dl);

	double * btbx = pfdaAlloc_d(*px**px**N, &dp);
	pfda_computebtb(btbx, N, Bx, M, px, nobs, dl);


	if(*kz){
		double * left = pfdaAlloc_d(*kz**kz,&dp);
		dsyrk_(&Upper, &Trans, kz, M, &dOne, z, M, &dzero, left, kz);
		dgemv_(&Trans, M, kz, &dOne, z, M, y, &one, &dzero, tz, &one);
		int sr =0, lwork = 10**kz;
		dsysv_(&Upper, kz, &one, left, kz, ip, tz, kz, dp, &lwork, &sr);
		if(sr)pfda_error("PFDA ERROR:dual_ca_i: computing tz dsysv returned:%d",sr);
	}

	pfda_s_i(tt, tf, gamma, Dg, gg, sigma, y,  nobs, M,  N,  kg,  Bt,  btbt, pt,  minV, dl, dp, ip);
	pfda_s_i(tx, tg, delta, Dd, dd, sigma, y,  nobs, M,  N,  kd,  Bx,  btbx, px,  minV, dl, dp, ip);

	for(int i=0;i<*N;i++)dger_(kg,kd,&dOne, gamma+i,N,delta+i,N, gd+i**kg**kd, kg);

	dual_bc_5(lambda, gg, gd, N, kg, kd, dl, dp, ip);

	double * Ry = pfdaAlloc_d(*M, &dp);
	dual_ca_resid( Ry, y, z, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, kz, kd, kg, pt, px, dl, dp);
	*sigma  = ddot_(M,Ry,&one,Ry,&one);
	*sigma /= *M;
}

/*!

@MEMORY
	 - dp = kg^2 + kd^2 + max(kg,kd)*10 + 2*kg*kd
	 - ip = max(kg,kd)
*/
void dual_ca_E1(
	double * const gamma,
	double * const delta,
	double * const gg,
	double * const gd,
	double * const dd,
	double * const Sgg,
	double * const Sgd,
	double * const Sdd,
	const double * const Ri,
	const double * const phi,
	const double * const psi,
	const double * const sigma,
	const double * const O1,
	const double * const O2,
	const double * const O3,
	const int * const ni,
	const int * const N,
	const int * const M,
	const int * const kg,
	const int * const kd,
	double * dp, int*ip)
{
	int kg2=*kg**kg, kgd=*kg**kd, kd2=*kd**kd;
	double si = 1/(*sigma);

	dcopy_(&kg2, O1,&one, Sgg, &one);
	dsyrk_(&Upper, &Trans, kg, ni, &si, phi, M, &dOne, Sgg, kg);

	dcopy_(&kgd, O2,&one, Sgd, &one);
	dgemm_(&Trans, &NoTrans, kg,kd,ni, &si, phi,M, psi,M, &dOne, Sgd, kg);

	dcopy_(&kd2, O3,&one, Sdd, &one);
	dsyrk_(&Upper, &Trans, kd, ni, &si, psi, M, &dOne, Sdd, kd);

	gen_symblock_solve(Sgg, Sgd,  Sdd, kg, kd, dp, ip);

	double * a  = pfdaAlloc_d(*kg,&dp);
	dgemv_(&Trans, ni, kg, &si, phi, M, Ri, &one, &dzero, a,&one);

	double * b  = pfdaAlloc_d(*kd,&dp);
	dgemv_(&Trans, ni, kd, &si, psi, M, Ri, &one, &dzero, b,&one);

	dsymv_(&Upper,   kg,    &dOne, Sgg, kg, a, &one, &dzero, gamma, N);
	dgemv_(&NoTrans, kg,kd, &dOne, Sgd, kg, b, &one, &dOne , gamma, N);

	dsymv_(&Upper,   kd,    &dOne, Sdd, kd, b, &one, &dzero, delta, N);
	dgemv_(&Trans,   kg,kd, &dOne, Sgd, kg, a, &one, &dOne , delta, N);

	dcopy_(&kg2,Sgg,&one,gg,&one);
	dsyr_(&Upper, kg, &dOne, gamma, N, gg, kg);
	pfda_fillsym(gg,kg, NULL);

	dcopy_(&kgd,Sgd,&one,gd,&one);
	dger_(kg,kd, &dOne, gamma, N, delta, N, gd, kg);

	dcopy_(&kd2,Sdd,&one,dd,&one);
	dsyr_(&Upper, kd, &dOne, delta, N, dd, kd);
	pfda_fillsym(dd,kd, NULL);
}

/*! Computes the E step estimating gamma, delta, Sgg,Sgd,Sdd,gg,gd,dd

@MEMORY
- dp = 2^kg^2 + 2*kd^2 + kg*kd+ max(kg,kd)*10 + 2*kg*kd+(1+kg+kd)*M
- ip = max(kg,kd)
*/
void dual_ca_E(
	double * const gamma,
	double * const delta,
	double * const gg,
	double * const gd,
	double * const dd,
	double * const Sgg,
	double * const Sgd,
	double * const Sdd,
	const double * const y,
	const double * const Z,
	const double * const Bt,
	const double * const Bx,
	const double * const tz,
	const double * const tt,
	const double * const tx,
	const double * const tf,
	const double * const tg,
	const double * const lambda,
	const double * const Dg,
	const double * const Dd,
	const double * const sigma,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const int * const dl,
	double * dp, int*ip)
{
	double * Mg = pfdaAlloc_d(*kg**kg,&dp); addDiag(kg, Dg, Mg, kg);
	double * Md = pfdaAlloc_d(*kd**kd,&dp); addDiag(kd, Dd, Md, kd);

	int kgd=*kg**kd;
	double * Mgd= pfdaAlloc_d(kgd,&dp);
	dgemm_(&NoTrans, &Trans, kg,kd,kg, &dOne, Mg, kg, lambda, kd, &dzero, Mgd, kg);

	gen_symblock_solve(Mg, Mgd,  Md, kg, kd, dp, ip);
	
	double * phi = pfdaAlloc_d(*M**kg,&dp);
	dgemm_(&NoTrans,&NoTrans, M, kg, pt, &dOne, Bt, M, tf, pt, &dzero, phi, M);

	double * psi = pfdaAlloc_d(*M**kd,&dp);
	dgemm_(&NoTrans,&NoTrans, M, kd, px, &dOne, Bx, M, tg, px, &dzero, psi, M);

	double * Ry = pfdaAlloc_d(*M,&dp);
	dual_ca_resid(Ry, y, Z, Bt, Bx, tz, tt, tx, NULL, NULL, NULL, NULL, nobs, N, M, kz, kg, kd, pt, px, dl, dp);
	for(int i=0,no=0;i<*N;i++){
		dual_ca_E1(
			gamma+i, delta+i,
			 gg+i**kg**kg,  gd+i**kg**kd,  dd+i**kd**kd,
			Sgg+i**kg**kg, Sgd+i**kg**kd, Sdd+i**kd**kd,
			Ry+no, phi+no,  psi+no,
			sigma,  Mg, Mgd, Md, nobs+i, N, M, kg, kd, dp, ip);
		no+=nobs[i];
	}
}

/*!  Computes the unpenalized coefficient vector

@MEMORY
- dp = M + kz^2 + 10*kz + M*max(kg, kd)
- ip = kz
*/
void dual_ca_unpenalized(
	double * const tz,
	const double * const y,
	const double * const Z,
	const double * const Bt,
	const double * const Bx,
	const double * const tt,
	const double * const tx,
	const double * const tf,
	const double * const tg,
	const double * const gamma,
	const double * const delta,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const int * const dl,
	double * dp, int*ip)
{
	double * Ry = pfdaAlloc_d(*M, &dp);
	dual_ca_resid( Ry, y, NULL, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, &zero, kd, kg, pt, px, dl, dp);
	
	if(Z && *kz){
		double * left = pfdaAlloc_d(*kz**kz,&dp);
		dsyrk_(&Upper, &Trans, kz, M, &dOne, Z, M, &dzero, left, kz);
		dgemv_(&Trans, M, kz, &dOne, Z, M, Ry, &one, &dzero, tz, &one);
		int sr =0, lwork = 10**kz;
		dsysv_(&Upper, kz, &one, left, kz, ip, tz, kz, dp, &lwork, &sr);
		if(sr)pfda_error("PFDA ERROR:dual_ca_i: computing tz dsysv returned:%d",sr);
	}
}

/*! Computes the estimates for the penalized parameter vectors (Not Principal Components)

@MEMORY
- dp =  10**p + p^2
-ip =  p
*/
void dual_ca_penalized_1(
	double * const tm,
	const double * const Ry,
	const double * const B,
	const double * const sigma,
	const double * const l,
	const double * const K,
	const int * const M,
	const int * const p,
	double * dp, int * ip)
{
	int p2=*p**p;
	double * left = pfdaAlloc_d(p2,&dp);
	double factor = *l**sigma;
	dsyrk_(&Upper, &Trans, p, M, &dOne, B, M, &dzero, left, p);
	daxpy_(&p2, &factor, K, &one, left, &one);
	
	dgemv_(&Trans, M, p, &dOne, B, M, Ry, &one, &dzero, tm, &one);
	
	int sr = 0,lwork = 10**p;
	dsysv_(&Upper, p,&one, left,  p, ip, tm, p, dp, &lwork, &sr);
}

/*! 

@MEMORY
- dp = 10*p + p^2 +M
- ip = p
*/
void dual_ca_penalized(
	double * const tt,
	double * const tx,
	const double * const y,
	const double * const Z,
	const double * const Bt,
	const double * const Bx,
	const double * const tz,
	const double * const tf,
	const double * const tg,
	const double * const gamma,
	const double * const delta,
	const double * const sigma,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const double * const lt,
	const double * const lx,
	const double * const Kt,
	const double * const Kx,
	const int * const dl,
	double * dp, int*ip)
{
	double * Ry = pfdaAlloc_d(*M,&dp);
	for(int i=0;i<*pt;i++)tt[i]=dzero;
	dual_ca_resid( Ry, y, Z, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, kz, kd, kg, pt, px, dl, dp);
	dual_ca_penalized_1( tt, Ry, Bt, sigma, lt, Kt, M, pt, dp, ip);
	
	for(int i=0;i<*px;i++)tx[i]=dzero;
	dual_ca_resid( Ry, y, Z, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, kz, kd, kg, pt, px, dl, dp);
	dual_ca_penalized_1( tx, Ry, Bx, sigma, lx, Kx, M, px, dp, ip);
}

/*!

@MEMORY 
 no memory requirements
*/
void dual_ca_princcomp_left(
	double * const left,
	const double * const Bt,
	const double * const sigma,
	const double * const gg,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kg,
	const int * const pt,
	const double * const lf,
	const double * const K,
	const int * const colnum)
{
	int p2 = *pt**pt, kg2 = *kg**kg;
	for(int i=0;i<p2;i++)left[i]=0;
	for(int i=0, no=0; i<*N;i++){
		dsyrk_(&Upper, &Trans, pt, nobs+i, gg+i*kg2+*colnum*(*kg+1), Bt+no, M, &dOne, left, pt);
		no+=nobs[i];
	}
	double factor = *lf**sigma;
	daxpy_(&p2,&factor, K, &one, left, &one);
}

/*! computes the right side of the equation for finding a column of Theta_f

@MEMORY 
- dp = M + N*(kg+kd) + M*max(kg,kd)
*/
void dual_ca_princcomp_right_tf(
	double * const right,
	double * const tf,
	const double * const Ry,
	const double * const Bt,
	const double * const Bx,
	const double * const tg,
	const double * const gamma,
	const double * const sigma,
	const double * const gg,
	const double * const gd,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const int * const colnum,
	double * dp)
{
	int kg2 = *kg**kg, kgd = *kg**kd;
	
	for(int i=0;i<*pt;i++)tf[i+*colnum**pt]=0.0;
	double * hg = pfdaAlloc_d(*N**kg,&dp);
	double * hd = pfdaAlloc_d(*N**kd,&dp);
	
	double * R2 = pfdaAlloc_d(*M,&dp); dcopy_(M,Ry, &one, R2,&one);
	for(int i=0, no=0; i<*N;i++){
		dcopy_(kg, gg+i*kg2+*colnum,kg, hg+i, N); 
		dcopy_(kd, gd+i*kgd+*colnum,kg, hd+i, N); 
		dscal_(nobs+i,gamma+i+*colnum**N, R2+no,&one);
		no+=nobs[i];
	}
	
	dual_ca_resid(R2, R2, NULL, Bt, Bx, NULL, NULL, NULL, tf, tg, hg, hd, nobs, N, M, NULL, kd, kg, pt, px, NULL, dp);
	dgemv_(&Trans, M, pt, &dOne, Bt, M, R2, &one, &dzero, right, &one);
}

/*! computes the right side of the equation for finding a column of Theta_g

it is iimportant to remember the distinction because here we are selection out columns not rows.
@MEMORY 
- dp = M + N*(kg+kd) + M*max(kg,kd)
*/
void dual_ca_princcomp_right_tg(
	double * const right,
	double * const tg,
	const double * const Ry,
	const double * const Bt,
	const double * const Bx,
	const double * const tf,
	const double * const delta,
	const double * const sigma,
	const double * const gd,
	const double * const dd,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const int * const colnum,
	double * dp)
{
	int kd2 = *kd**kd, kgd = *kg**kd;
	
	for(int i=0;i<*px;i++)tg[i+*colnum**px]=0.0;
	// these represent nearly the same thins as the tf version but are computed reverse
	double * hg = pfdaAlloc_d(*N**kg,&dp);
	double * hd = pfdaAlloc_d(*N**kd,&dp);
	
	double * R2 = pfdaAlloc_d(*M,&dp); dcopy_(M,Ry, &one, R2,&one);
	for(int i=0, no=0; i<*N;i++){
		dcopy_(kd, dd+i*kd2+*colnum**kd,&one, hd+i, N); 
		dcopy_(kg, gd+i*kgd+*colnum**kg,&one, hg+i, N); 
		dscal_(nobs+i,delta+i+*colnum**N, R2+no,&one);
		no+=nobs[i];
	}
	
	dual_ca_resid(R2, R2, NULL, Bt, Bx, NULL, NULL, NULL, tf, tg, hg, hd, nobs, N, M, NULL, kd, kg, pt, px, NULL, dp);
	dgemv_(&Trans, M, px, &dOne, Bx, M, R2, &one, &dzero, right, &one);
}

/*! computes the new estimates of Theta_f and Theta_g

@MEMORY
 - dp = M + p^2 + p + max( 10*pt ,  M + N*(kg+kd) + M*max(kg,kd) )
 - ip = p
*/
void dual_ca_princcomp(
	double * const tf,
	double * const tg,
	const double * const y,
	const double * const Z,
	const double * const Bt,
	const double * const Bx,
	const double * const tz,
	const double * const tt,
	const double * const tx,
	const double * const gamma,
	const double * const delta,
	const double * const sigma,
	const double * const gg,
	const double * const gd,
	const double * const dd,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const double * const lf,
	const double * const lg,
	const double * const Kt,
	const double * const Kx,
	const int * const dl,
	double * dp, int*ip)
{
	double pp =MMAX(*pt,*px);double * p = &pp;
	double * Ry = pfdaAlloc_d(*M,&dp);
	dual_ca_resid( Ry, y, Z, Bt, Bx, tz, tt, tx, NULL, NULL, gamma, delta, nobs, N, M, kz, kd, kg, pt, px, NULL, dp);
	
	double * left = pfdaAlloc_d(*p**p,&dp);
	double * right = pfdaAlloc_d(*p,&dp);
	for(int j=0; j<*kg;j++){
		dual_ca_princcomp_left( left, Bt, sigma, gg, nobs, N, M, kg, pt, lf, Kt, &j);
		
		dual_ca_princcomp_right_tf( right, tf, Ry, Bt, Bx, tg, gamma, sigma, gg, gd, nobs, N, M, kg, kd, pt, px, &j, dp);
		
		int lwork = 10**p,sr=0;
		dsysv_(&Upper, pt, &one, left, pt, ip, right, pt, dp, &lwork, &sr);
		u_orthogonalize(right,  tf, pt, &j);
		dcopy_(pt,right,&one,tf+j**pt,&one);
	}
	for(int j=0; j<*kd;j++){
		dual_ca_princcomp_left( left, Bx, sigma, dd, nobs, N, M, kd, px, lg, Kx, &j);
		
		dual_ca_princcomp_right_tg( right, tg, Ry, Bt, Bx, tf, delta, sigma, gd, dd, nobs, N, M, kg, kd, pt, px, &j, dp);
		
		int lwork = 10**p,sr=0;
		dsysv_(&Upper, px, &one, left, px, ip, right, px, dp, &lwork, &sr);
		u_orthogonalize(right,  tg, px, &j);
		dcopy_(px,right,&one,tg+j**px,&one);
	}
}

/*! estimates the variance parameters of the model

also updates variance associated parameters gamma and delta.
@MEMORY
- dp = 2*kg^2 + 2*kd^2 + M + 10*k + 9*p + 3*p^2 + N*k + M*k
- ip = 6*p 
*/
void dual_ca_variances(
	      double * const tf,
	      double * const tg,
	      double * const gamma,
	      double * const delta,
	      double * const lambda,
	      double * const Dg,
	      double * const Dd,
	      double * const sigma,
	const double * const y,
	const double * const Z,
	const double * const Bt,
	const double * const Bx,
	const double * const tz,
	const double * const tt,
	const double * const tx,
	const double * const gg,
	const double * const gd,
	const double * const dd,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const int * const dl,
	double * dp, int* ip)
{
	double * Ry = pfdaAlloc_d(*M, &dp);
	dual_ca_resid( Ry, y, Z, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, kz, kd, kg, pt, px, dl, dp);
	*sigma  = ddot_(M,Ry,&one,Ry,&one);
	*sigma /= *M;
	
	dual_bc_5(lambda, gg, gd, N, kg, kd, dl, dp, ip);

	int kg2=*kg**kg, kd2 = *kd**kd;
	double Ninv = 1/(double)*N;
	double * sumgg = pfdaAlloc_d(kg2,&dp);
	double * sumdd = pfdaAlloc_d(kd2,&dp);
	for(int i=0;i<*N;i++){
		daxpy_(&kg2,&Ninv, gg+i*kg2,&one,sumgg,&one);
		daxpy_(&kd2,&Ninv, dd+i*kd2,&one,sumdd,&one);
	}	
	
	double * transa = pfdaAlloc_d(kg2,&dp);
	double * transb = pfdaAlloc_d(kd2,&dp);
	gen_orthog( tf, gamma, Dg, transa, sumgg, N, kg, pt, dl, dp, ip);
	gen_orthog( tg, delta, Dd, transb, sumdd, N, kd, px, dl, dp, ip);
	
	pfdaDual_m5_2(lambda, transa, kg, transb, kd, dl,dp,ip);
}

/*! Contains the primary loop for the EM algorithm for the calicum model.

The calcium model is a nick name for a single response with un penalized predictors as well as two sets of principle components from two predeictor variables.
@MEMORY
- dp = kz + kg + kd + pt +px + pt*kg + px*kd + kg*kd + max of the steps
- ip = max of steps = max(6*p, kz)
*/
void dual_ca_core(
	double * const y,
	double * const Z,
	double * const Bt,
	double * const Bx,
	double * const tz,
	double * const tt,
	double * const tx,
	double * const tf,
	double * const tg,
	double * const gamma,
	double * const delta,
	double * const lambda,
	double * const Dg,
	double * const Dd,
	double * const sigma,
	double * const gg,
	double * const gd,
	double * const dd,
	double * const Sgg,
	double * const Sgd,
	double * const Sdd,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kz,
	const int * const kg,
	const int * const kd,
	const int * const pt,
	const int * const px,
	const double * const lt,
	const double * const lx,
	const double * const lf,
	const double * const lg,
	const double * const Kt,
	const double * const Kx,
	const double * const minV,
	      int * const maxI,
	      double * const tol,
	const int * const dl, double * dp, int * ip)
{
	pfda_debug_cdl(debugnum_dual_ca_inputs){
		pfda_debug_arg(*lt);
		pfda_debug_arg(*lx);
		pfda_debug_arg(*lf);
		pfda_debug_arg(*lg);
		pfda_debug_arg(Kt[1]);
		pfda_debug_arg(Kx[1]);
		pfda_debug_arg(*minV);
		pfda_debug_arg(*maxI);
		pfda_debug_arg(*tol);
		pfda_debug_arg(*dl);
	}
	dual_ca_i( tz, tt, tx, tf, tg, gamma, delta, lambda, Dg, Dd, sigma, gg, gd, dd, y, Z, Bt, Bx, nobs, N, M, kz, kg, kd, pt, px, minV, dl, dp, ip);
	int I=0, kgd = *kg**kd, pg = *pt**kg, pd=*px**kd;
	double cc=0;
	double 	sigma_old=0;	
	double * tz_old     =pfdaAlloc_d(*kz,     &dp);
	double * tt_old     =pfdaAlloc_d(*pt,     &dp);
	double * tx_old     =pfdaAlloc_d(*px,     &dp);
	double * tf_old     =pfdaAlloc_d(*pt**kg, &dp);
	double * tg_old     =pfdaAlloc_d(*px**kd, &dp);
	double * lambda_old =pfdaAlloc_d(*kg**kd, &dp);
	double * Dg_old     =pfdaAlloc_d(*kg,     &dp);
	double * Dd_old     =pfdaAlloc_d(*kd,     &dp);
	while(I++<*maxI)/*limited loop with exit*/{
		{ //setup for convergence
			sigma_old = *sigma;
			dcopy_(kz  , tz    , &one, tz_old    , &one);
			dcopy_(pt  , tt    , &one, tt_old    , &one);
			dcopy_(px  , tx    , &one, tx_old    , &one);
			dcopy_(&pg , tf    , &one, tf_old    , &one);
			dcopy_(&pd , tg    , &one, tg_old    , &one);
			dcopy_(&kgd, lambda, &one, lambda_old, &one);
			dcopy_(kg  , Dg    , &one, Dg_old    , &one);
			dcopy_(kd  , Dd    , &one, Dd_old    , &one);
		}
		dual_ca_E( gamma, delta, gg, gd, dd, Sgg, Sgd, Sdd, y, Z, Bt, Bx, tz, tt, tx, tf, tg, lambda, Dg, Dd, sigma, nobs, N, M, kz, kg, kd, pt, px, dl, dp, ip);
		dual_ca_unpenalized( tz, y, Z, Bt, Bx, tt, tx, tf, tg, gamma, delta,  nobs, N, M, kz, kg, kd, pt, px, dl, dp, ip);
		dual_ca_penalized( tt, tx, y, Z, Bt, Bx, tz, tf, tg, gamma, delta, sigma, nobs, N, M, kz, kg, kd, pt, px, lt, lx, Kt, Kx, dl, dp, ip);
		dual_ca_princcomp( tf, tg, y, Z, Bt, Bx, tz, tt, tx, gamma, delta, sigma, gg, gd, dd, nobs, N, M, kz, kg, kd, pt, px, lf, lg, Kt, Kx, dl, dp, ip);
		dual_ca_variances( tf, tg, gamma, delta, lambda, Dg, Dd, sigma, y, Z, Bt, Bx, tz, tt, tx, gg, gd, dd, nobs, N, M, kz, kg, kd, pt, px, dl, dp, ip);
		{ // Compute Convergence Criteria 
			cc = fabs(sigma_old-*sigma)/MMAX(*sigma,*minV);

			double cctz=0;
			if(tz && kz && *kz){
				daxpy_(kz, &mOne, tz, &one, tz_old, &one);
				for(int i=0;i<*kz;i++){
					cctz+=fabs(tz_old[i])/MMAX(fabs(tz[i]),*minV);
				}
				cc += cctz/ (double)*kz;
			}
			
			daxpy_(pt, &mOne, tt, &one, tt_old, &one);
			double cctt=0;
			for(int i=0;i<*pt;i++){
				cctt+=fabs(tt_old[i])/MMAX(fabs(tt[i]),*minV);
			}
			cc+=cctt/(double)(*pt);

			daxpy_(px, &mOne, tx, &one, tx_old, &one);
			double cctx=0;
			for(int i=0;i<*px;i++){
				cctx+=fabs(tx_old[i])/MMAX(fabs(tx[i]),*minV);
			}
			cc+=cctx/(*px);

			double cctf=0;
			daxpy_(&pg, &mOne, tf, &one, tf_old, &one);
			for(int i=0;i<pg;i++)cctf+=fabs(tf_old[i])/MMAX(fabs(tf[i]),*minV);
			cc+=cctf/(double)pg;

			double cctg=0;
			daxpy_(&pd, &mOne, tg, &one, tg_old, &one);
			for(int i=0;i<pd;i++)cctg+=fabs(tg_old[i])/(MMAX(fabs(tg[i]),*minV));
			cc+=cctg/(double)pd;

			daxpy_(&kgd, &mOne, lambda, &one, lambda_old, &one);
			double ccL=0;
			for(int i=0;i<kgd;i++)ccL+=fabs(lambda_old[i])/MMAX(fabs(lambda[i]),*minV);
			cc+=ccL/(double)kgd;
			
			daxpy_(kg, &mOne, Dg, &one, Dg_old, &one);
			double ccDg=0;
			for(int i=0;i<*kg;i++)ccDg+=fabs(Dg_old[i])/MMAX(fabs(Dg[i]),*minV);
			cc+=ccDg/ (double)*kg;

			daxpy_(kd, &mOne, Dd, &one, Dd_old, &one);
			double ccDd=0;
			for(int i=0;i<*kd;i++)ccDd+=fabs(Dd_old[i])/MMAX(fabs(Dd[i]),*minV);
			cc+=ccDd/ (double)*kd;

			cc/=9;
			///* Debug */
			pfda_debug_cdl(debugnum_dual_criteria_components){
				pfda_debug_msg("Convergence criteria components:\n");
				pfda_debug_msg("sigma:  \t%5.5g\n", fabs(sigma_old-*sigma));
				pfda_debug_msg("tz:   \t%5.5g\n", cctz);
				pfda_debug_msg("tt:   \t%5.5g\n", cctt);
				pfda_debug_msg("tx:   \t%5.5g\n", cctx);
				pfda_debug_msg("Tf:   \t%5.5g\n", cctf);
				pfda_debug_msg("Tg:   \t%5.5g\n", cctg);
				pfda_debug_msg("lambda:    \t%5.5g\n", ccL);
				pfda_debug_msg("Dg:   \t%5.5g\n", ccDg);
				pfda_debug_msg("Dd:   \t%5.5g\n", ccDd);
				fflush(stdout);
			}
			pfda_debug_cdl(debugnum_dual_criteria){pfda_debug_arg(cc);fflush(stdout);}
			if(cc < *tol){ //break when acheiving convergence
				pfda_debug_cdl(debugnum_dual){pfda_debug_msg("Criteria met leaving loop.%g\n");fflush(stdout);}
				break;
			}	
		}
	}
	if(*maxI<=I){pfda_error("EM-algorithm did not converge");}
	dual_ca_resid(y, y, Z, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, kz, kd, kg, pt, px, dl, dp);
	*tol=cc;
	*maxI = I;	
}


