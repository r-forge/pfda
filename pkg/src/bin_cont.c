#include "pfda_single.h"
#include "pfda_gen.h"
#include "pfda_dual.h"
#include "pfda_const.h"
#include "pfda_memory.h"
#include "pfda_utils.h"
#include "pfda_error.h"
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

/*! Find the initial values for the half binary half continuous case

\group dual binary/continuous
@MEMORY
	- dp = 2*p^2 + 2*M + M*k+ p*N + 8*p
	- ip = 6*p 
*/
void dual_bc_i(
	      double * const tm,
	      double * const tn,
	      double * const tf,
	      double * const tg,
	      double * const alpha,
	      double * const beta,
	      double * const lambda,
	      double * const Da,
	      double * const Db,
	      double * const s_xi,
	      double * const aa,
	      double * const bb,
	const int    * const y,
	const double * const z,
	const double * const B,
	const int    * const nobs,
	const double * const btb,
	const int    * const N,
	const int    * const M,
	const int    * const ka,
	const int    * const kb,
	const int    * const p,
	const double * const minV,
	const int * const dl,
	double * dp, int*ip)
{
	double * w = pfdaAlloc_d(*M,&dp);
	for(int i=0;i<*M;i++)w[i]=(double)y[i];
	double s_eps=0;
	pfda_s_i(tm, tf, alpha, Da, aa, &s_eps, w, nobs, M, N, ka, B, btb, p,  minV, dl, dp, ip);
	pfda_s_i(tn, tg, beta , Db, bb, s_xi, z, nobs, M, N, kb, B, btb, p,  minV, dl, dp, ip);
	int ka2 = *ka**ka;
	double * sum_aa = pfdaAlloc_d(ka2,&dp);
	double * sum_ab  = lambda;
	
	dsyrk_(&Upper,&Trans,ka, N, &dOne, alpha, N, &dzero, sum_aa, ka);
	dgemm_(&Trans, &NoTrans, ka, kb, N, &dOne, alpha, N, beta, N, &dzero, sum_ab, ka);
	
	int svr=0;
	int * ipiv = pfdaAlloc_i(*N,&ip);
	int lwork = 8**N;
	dsysv_(&Upper, ka, kb, sum_aa, ka, ipiv, sum_ab, ka, dp, &lwork,&svr);
	pfda_transpose(lambda,*ka, *kb,dl,dp);
}

/*!  computes the  following \f$ \Sigma_{\alpha\alpha}B_i^T\Theta_f^T Ry_i+\Sigma_{\alpha\beta}B_i^T\Theta_f Rz_i/\sigma_\xi\f$

@MEMORY
-dp = ka +kb
*/
void dual_bc_1a(
	      double * const mu,
	const double * const Ry,
	const double * const Rz,
	const double * const phi,
	const double * const psi,
	const double * const sxi,
	const double * const Sa,
	const double * const Sab,
	const int    * const M,
	const int    * const ni,
	const int    * const ka,
	const int    * const kb,
	const int * const dl,
	double*dp
){
	double * tmpa=pfdaAlloc_d(*ka, &dp);
	double * tmpb=pfdaAlloc_d(*kb, &dp);
	dgemv_(&Trans, ni,ka,&dOne, phi, M, Ry, &one, &dzero, tmpa, &one);
	dgemv_(&Trans, ni,kb,&dOne, psi, M, Rz, &one, &dzero, tmpb, &one);
	
	dsymv_(&Upper, ka, &dOne, Sa, ka, tmpa, &one, &dzero, mu,&one);  
	
	double sxi_inv = 1/(*sxi);
	dgemv_(&NoTrans, ka, kb, &sxi_inv, Sab, ka, tmpb, &one, &dOne, mu,&one);
}

/*!  computes the  following \f$ \Sigma_{\alpha\beta}^T B_i^T\Theta_f^T Ry_i+\Sigma_{\beta\beta}B_i^T\Theta_f Rz_i /\sigma_\xi\f$

@MEMORY
-dp = ka +kb
*/
void dual_bc_1b(
	      double * const mu,
	const double * const Ry,
	const double * const Rz,
	const double * const phi,
	const double * const psi,
	const double * const sxi,
	const double * const Sab,
	const double * const Sbb,
	const int    * const M,
	const int    * const ni,
	const int    * const ka,
	const int    * const kb,
	const int * const dl,
	double*dp
){
	double * tmpa=pfdaAlloc_d(*ka, &dp);
	double * tmpb=pfdaAlloc_d(*kb, &dp);
	dgemv_(&Trans, ni,ka,&dOne, phi, M, Ry, &one, &dzero, tmpa, &one);
	dgemv_(&Trans, ni,kb,&dOne, psi, M, Rz, &one, &dzero, tmpb, &one);
	
	double sxi_inv = 1/(*sxi);
	dsymv_(&Upper, kb, &sxi_inv, Sbb, kb, tmpb, &one, &dzero, mu,&one);  
	
	dgemv_(&Trans, ka, kb, &dOne, Sab, ka, tmpa, &one, &dOne, mu,&one);
}


/*!   Computes aa, ab, bb
 

@MEMORY
-dp = 3ka +3kb + ka*ni + kb*ni + ni*max(ka,kb)
*/
void dual_bc_1cde(
	      double * const aa,
	      double * const ab,
	      double * const bb,
	const double * const Rz,
	const double * const wi,
	const double * const rho,
	const double * const phi,
	const double * const psi,
	const double * const wwi,
	const double * const tf,
	const double * const tg,
	const double * const sxi,
	const double * const Saa,
	const double * const Sab,
	const double * const Sbb,
	const int    * const M,
	const int    * const ni,
	const int    * const ka,
	const int    * const kb,
	const int * const dl,
	double*dp)
{
	double * S1 = pfdaAlloc_d(*ka, &dp);
	dual_bc_1a(S1,rho,Rz,phi,psi,sxi,Saa,Sab,M,ni,ka,kb,dl,dp);

	double * S2 = pfdaAlloc_d(*kb, &dp);
	dual_bc_1b(S2,rho,Rz,phi,psi,sxi,Sab,Sbb,M,ni,ka,kb,dl,dp);

	double * S3 = pfdaAlloc_d(*ka**ni, &dp);  //defined as transpose of what is foudn in R version
	dsymm_(&Right, &Upper, ni,ka, &dOne, Saa, ka, phi, M, &dzero, S3, ni);

	double * S4 = pfdaAlloc_d(*ka, &dp);
	dgemv_(&Trans, ni, ka, &dOne, S3, ni, wi, &one, &dzero, S4, &one);
	
	double * S5 = pfdaAlloc_d(*ni**kb, &dp);
	dgemm_(&Trans, &Trans, kb, ni, ka, &dOne, Sab, ka, phi, M, &dzero, S5, kb);
	
	double * S6 = pfdaAlloc_d(*kb, &dp);
	dgemv_(&NoTrans, kb, ni, &dOne, S5, kb, wi, &one, &dzero, S6, &one);
	
	int ka2=*ka**ka;
	pfda_matrix_inner_quadratic_form( aa,S3,ka,ni,wwi,ni,dl,dp);
	dsyr2_(&Upper, ka, &dOne, S1, &one, S4, &one, aa, ka);
	dsyr_(&Upper, ka, &dOne, S1, &one, aa, ka);
	daxpy_(&ka2,&dOne,Saa,&one,aa,&one);

	int kb2 = *kb**kb;
	pfda_matrix_outer_quadratic_form( bb,S5,kb,kb,wwi,ni,dl,dp);
	dsyr2_(&Upper, kb, &dOne, S2, &one, S6, &one, bb, kb);
	dsyr_(&Upper, kb, &dOne, S2, &one, bb, kb);
	daxpy_(&kb2,&dOne,Sbb,&one,bb,&one);
	
	int kab = *ka**kb;
	double * tmp=pfdaAlloc_d(*ni**ka,&dp);
	dsymm_(&Left, &Upper, ni,ka, &dOne, wwi, ni, S3, ni,&dzero, tmp, ni);
	dgemm_(&Trans, &Trans, ka, kb, ni, &dOne, tmp, ni, S5, kb, &dzero, ab, ka);
	dger_(ka, kb, &dOne,S1, &one, S6, &one, ab, ka);
	dger_(ka, kb, &dOne,S4, &one, S2, &one, ab, ka);
	dger_(ka, kb, &dOne,S1, &one, S2, &one, ab, ka);
	daxpy_(&kab,&dOne, Sab, &one,ab,&one);
}


/*!   E step for binary/continuous model

@MEMORY
- dp = M*(ka + kb + 3) + max(7 * max(ka,kb)^2 , 3ka +3kb + (ka + kb+ max(ka,kb))*ni)
 - ip = max(ka,kb)
*/
void dual_bc_1(
	      double * const alpha,
	      double * const beta,
	      double * const aa,
	      double * const ab,
	      double * const bb,
	      double * const Saa,
	      double * const Sab,
	      double * const Sbb,
	const double * const z,
	const double * const B,
	const double * const w,
	const double * const ww,
	const double * const tm,
	const double * const tn,
	const double * const tf,
	const double * const tg,
	const double * const lambda,
	const double * const Da,
	const double * const Db,
	const double * const sxi,
	const int    * const nobs,
	const int    * const N,
	const int    * const M,
	const int    * const ka,
	const int    * const kb,
	const int    * const p,
	const int * const dl,
	double*dp, int * ip)
{
	double * phi = pfdaAlloc_d(*M**ka,&dp);
	dgemm_(&NoTrans,&NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);
	
	double * psi = pfdaAlloc_d(*M**kb,&dp);
	dgemm_(&NoTrans,&NoTrans, M, kb, p, &dOne, B, M, tg, p, &dzero, psi, M);
	
	double * rho = pfdaAlloc_d(*M, &dp);
	dgemv_(&NoTrans, M,p, &mOne, B, M, tm, &one, &dzero, rho, &one);
	
	double * Rw = pfdaAlloc_d(*M, &dp);
	dcopy_(M,w,&one,Rw,&one);
	daxpy_(M,&dOne,rho,&one,Rw,&one);
	
	double * Rz = pfdaAlloc_d(*M, &dp);
	dcopy_(M,z,&one,Rz,&one);
	dgemv_(&NoTrans, M,p, &mOne, B, M, tn, &one, &dOne, Rz, &one);
	
	double * a = pfdaAlloc_d(*ka,&dp);
	double * b = pfdaAlloc_d(*kb,&dp);
	int no = 0, nno=0;
	for(int i=0;i<*N;i++){
		dual_gen_sigmas( Saa+i**ka**ka, Sab+i**ka**kb, Sbb+i**kb**kb, phi+no, psi+no, lambda, Da, Db, &dOne, sxi, M, nobs+i, ka, kb, dl, dp, ip);
		dual_bc_1a( a, Rw+no, Rz+no, phi+no, psi+no, sxi, Saa+i**ka**ka , Sab+i**ka**kb, M, nobs+i, ka, kb, dl, dp);
		dcopy_(ka,a,&one,alpha+i,N);
		dual_bc_1b( b, Rw+no, Rz+no, phi+no, psi+no, sxi, Sab+i**ka**kb, Sbb+i**kb**kb, M, nobs+i, ka, kb, dl, dp);
		dcopy_(kb,b,&one,beta+i,N);
		dual_bc_1cde(aa+i**ka**ka, ab+i**ka**kb, bb+i**kb**kb, Rz+no, w+no, rho+no, phi+no, psi+no, ww+nno, tf, tg, sxi, Saa+i**ka**ka, Sab+i**ka**kb, Sbb+i**kb**kb, M, nobs+i, ka, kb, dl, dp);
		no  += nobs[i];
		nno += nobs[i]*nobs[i];
	}
}

/*! computes the estimate of 
@Memory
-dp = M + M*kb + 2*kb^2
*/
void dual_bc_2(
	      double * const sxi,
	const double * const z,
	const double * const B,
	const double * const tn,
	const double * const tg,
	const double * const beta,
	const double * const Sbb,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const kb,
	const int * const p,
	const int * const dl,
	double * dp)
{
	pfda_debug_msg("z[1]=%g\n",*z);
	pfda_debug_msg("B[1]=%g\n",*B);
	pfda_debug_msg("tn[1]=%g\n",*tn);
	pfda_debug_msg("tg[1]=%g\n",*tg);
	pfda_debug_msg("beta[1]=%g\n",*beta);
	pfda_debug_msg("Sbb[1]=%g\n",*Sbb);
	pfda_debug_msg("nobs[1]=%d\n",*nobs);
	pfda_debug_msg("N[1]=%d\n",*N);
	pfda_debug_msg("M[1]=%d\n",*M);
	pfda_debug_msg("kb[1]=%d\n",*kb);
	pfda_debug_msg("p[1]=%d\n",*p);
	pfda_debug_msg("dl[1]=%d\n",*dl);
	
	double minimum_variance = 1e-4;
	pfda_m1(sxi,z,nobs,M,N,kb,B,p,&minimum_variance,tn,tg,beta,Sbb,dl, dp);
}

/*! computing tm, & tn
@MEMORY
	- dp =  p^2+ M+ M*max(ka,kb)
*/
void dual_bc_3(
	      double * const tm,
	      double * const tn,
	const double * const z,
	const double * const B,
	const double * const w,
	const double * const tf,
	const double * const tg,
	const double * const alpha,
	const double * const beta,
	const double * const sxi,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const ka,
	const int * const kb,
	const int * const p,
	const double * const lm,
	const double * const ln,
	const double * const K,
	const int * const dl,
	double * dp)
{
	pfda_m2(tm, w, nobs,M, N, ka, B, p, lm, K, tf, &dOne, alpha, dl, dp);
	pfda_m2(tn, z, nobs,M, N, kb, B, p, ln, K, tg, sxi  , beta , dl, dp);
}

/*! computes tf and tg
@MEMORY
	- dp = M + max(ka,kb)^2 * N + p^2 + p 
	- ip = p
*/
void dual_bc_4(
	      double * const tf,
	      double * const tg,
	const double * const z,
	const double * const B,
	const double * const w,
	const double * const tm,
	const double * const tn,
	const double * const alpha,
	const double * const beta,
	const double * const sxi,
	const double * const aa,
	const double * const bb,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const ka,
	const int * const kb,
	const int * const p,
	const double * const lf,
	const double * const lg,
	const double * const K,
	const double * const btb,
	const int * const dl,
	double * dp, int* ip
)
{
	pfda_m3_core(tf, aa, w, nobs, M, N, ka, B, p, lf, K, tm, &dOne, alpha, btb, dl, dp, ip);
	pfda_m3_core(tg, bb, z, nobs, M, N, kb, B, p, lg, K, tn, sxi  , beta , btb, dl, dp, ip);
}

/*! Computes  estimate of lambda

@MEMORY
	- dp =  ka^2 + ka*kb + 10*max(ka,kb)
	- ip = ka
*/
void dual_bc_5(
	double * lambda,
	const double * const aa,
	const double * const ab,
	const int * const N,
	const int * const ka,
	const int * const kb,
	const int * const dl,double * dp, int * ip)
{
	int ka2=*ka**ka, kab = *ka**kb;
	double * sumaa = pfdaAlloc_d(ka2,&dp);
	double * sumab = lambda; //ALIASING
	for(int i=0;i<kab;i++)sumab[i]=0;
	for(int i=0;i<*N;i++){
		daxpy_(&ka2,&dOne, aa+i*ka2,&one,sumaa,&one);
		daxpy_(&kab,&dOne, ab+i*kab,&one,sumab,&one);
	}
	int sr=0;
	int lwork = 10**MMAX(ka,kb);
	dsysv_(&Upper, ka, kb, sumaa, ka,ip,sumab,ka,dp, &lwork,&sr);
	pfda_transpose(lambda, *ka, *kb,dl, dp);
}

/*! Computes Da,and Db, updates, tf,tg,alpha,beta, & lambda,

@MEMORY
	- dp = 4*p^2 + 2*p +  max(  ka, kb, 8 ) * p + ka*kb
	- ip = 6*p
*/
void dual_bc_6(
	      double * const tf,
	      double * const tg,
	      double * const alpha,
	      double * const beta,
	      double * const lambda,
	      double * const Da,
	      double * const Db,
	const double * const aa,
	const double * const bb,
	const int * const N,
	const int * const ka,
	const int * const kb,
	const int * const p,
	const int * const dl,
	double * dp, int* ip)
{
	int ka2=*ka**ka, kb2 = *kb**kb;
	double Ninv = 1/(double)*N;
	double * sumaa = pfdaAlloc_d(ka2,&dp);
	double * sumbb = pfdaAlloc_d(kb2,&dp);
	for(int i=0;i<*N;i++){
		daxpy_(&ka2,&Ninv, aa+i*ka2,&one,sumaa,&one);
		daxpy_(&kb2,&Ninv, bb+i*kb2,&one,sumbb,&one);
	}	
	
	double * transa = pfdaAlloc_d(ka2,&dp);
	double * transb = pfdaAlloc_d(kb2,&dp);
	gen_orthog( tf, alpha, Da, transa, sumaa, N, ka, p, dl, dp, ip);
	gen_orthog( tg, beta , Db, transb, sumbb, N, kb, p, dl, dp, ip);
	
	pfdaDual_m5_2(lambda, transa, ka, transb, kb, dl,dp,ip);
}

/*! simulates W|Y,Z

@MEMORY
	- dp = ka + ka^2 + ka*kb +kb^2 + (7 * max(ka,kb)^2) + ka + (ka+kb) + ka + 1 + ka + 1
		  <= 2 + 5*ka + kb + 10 * k^2
	- ip = max(ka,kb)
*/
void dual_bc_genw(
	double       * const w_sim,
	int    const * const y,
	double const * const Rz,
	double const * const Rw,
	double const * const pi,
	double       * const phi,  // must be editable, but changes are reverted after.
	double const * const psi,
	double const * const lambda,
	double const * const Da,
	double const * const Db,
	int    const * const ni,
	int    const * const M, 
	int    const * const ka,
	int    const * const kb,
	int    const * const kr,
	int    const * const p, 
	int    const * const j,
	int const * const dl,
	double * dp, int * ip)
{
	double * old_phi_row = pfdaAlloc_d(*ka,&dp);
	dcopy_(ka, phi+*j, M, old_phi_row, &one);
	for(int i=0;i<*ka;i++)phi[*j+i**M]=dzero;
	
	double * Saa = pfdaAlloc_d(*ka**ka,&dp);
	double * Sab = pfdaAlloc_d(*ka**kb,&dp);
	double * Sbb = pfdaAlloc_d(*kb**kb,&dp);
	dual_gen_sigmas(Saa, Sab, Sbb, phi, psi, lambda, Da, Db, &dOne, &dOne, M, ni, ka, kb, dl, dp, ip);
	
	double * mu = pfdaAlloc_d(*ka,&dp);
	dual_bc_1a(mu, Rw, Rz, phi, psi, &dOne, Saa, Sab, M, ni, ka, kb, dl, dp);
	
	double * s = pfdaAlloc_d(one, &dp);
	pfda_matrix_inner_quadratic_form(s, old_phi_row, &one, ka, Saa, ka, dl, dp);
	*s += dOne;
	
	double * a = pfdaAlloc_d(one, &dp);
	*a = ddot_(ka, old_phi_row, &one, mu, &one);
	*a += pi[*j];
	
	if(y[*j]){
		double c = -*a/(*s);
		pfda_gen_truncnorm( w_sim, kr, &c, dl);
		for(int k=0;k<*kr;k++){
			w_sim[k] *= *s;
			w_sim[k] += *a;
		}
	} else {
		double c = *a/(*s);
		pfda_gen_truncnorm( w_sim, kr, &c, dl);
		for(int k=0;k<*kr;k++){
			w_sim[k] *= -*s;
			w_sim[k] += *a;
		}
	}

	dcopy_(ka, old_phi_row, &one, phi+*j, M);
}

/*!	computes the w and ww for a subject

\ingroup binary_single
@MEMORY:
	- dp length 	=  ni*kr + kr +   (2 + 5*ka + kb + 10 * k^2)
	- ip length = k
*/
void dual_bc_w_1(
	double * const w,
	double * const ww,
	int    const * const y,
	double const * const Rw,
	double const * const Rz,
	double const * const pi,
	double       * const phi,
	double const * const psi,
	double const * const lambda,
	double const * const Da,
	double const * const Db,
	int    const * const ni,
	int    const * const M,
	int    const * const ka,
	int    const * const kb,
	double const * const weight,
	int    const * const kr,
	int    const * const p,
	int const * const dl, double * dp, int * ip)
{
/* 
	double * Rw = pfdaAlloc_d(*ni,&dp);
	pfda_computeResid(Rw,w,ni, M, NULL, ka, B, p, tm, NULL, NULL, dl, dp);
	
	double * Rz = pfdaAlloc_d(*ni,&dp);
	pfda_computeResid(Rz,z,ni, M, NULL, kb, B, p, tn, NULL, NULL, dl, dp);

	double * phi = pfdaAlloc_d(*M**ka,&dp);
	dgemm_(&NoTrans,&NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);
	
	double * psi = pfdaAlloc_d(*M**kb,&dp);
	dgemm_(&NoTrans,&NoTrans, M, kb, p, &dOne, B, M, tg, p, &dzero, psi, M);
	
	double * pi = pfdaAlloc_d(*M, &dp);
	dgemv_(&NoTrans, M,p, &dOne, B, M, tm, &one, &dzero, pi, &one);
 */	
	double * wsim = pfdaAlloc_d(*ni**kr, &dp);
	for(int j=0;j<*ni;j++)
		dual_bc_genw(wsim+j**kr, y, Rz, Rw, pi, phi, psi, lambda, Da, Db, ni, M,  ka, kb, kr, p,  &j, dl, dp, ip);
	
	double krinv = 1.0/(double)*kr;
	double * V1 = pfdaAlloc_d(*kr,&dp);
	for(int i=0;i<*kr;i++)V1[i]=krinv;
	
	double cweight = (1-*weight);
	dgemv_(&Trans, kr, ni, weight, wsim, kr, V1, &one, &cweight, w, &one);

	double w2 = *weight/(double)*kr;
	dsyrk_(&Upper, &Trans, ni, kr, &w2, wsim, kr, &cweight, ww, ni);
	pfda_fillsym (ww,ni,dl);
}

/*! estimates by simulation \f[ E(W|Y,Z) \f] and \f[ E(WW|Y,Z) \f] for all subjects.
@ MEMORY
	- dp = M*(3+ka+kb) +   (ni*kr + kr + 2 + 5*ka + kb + 10 * k^2)
	- ip = k
*/
void dual_bc_w(
	double * const w,
	double * const ww,
	int    const * const y,
	double const * const z,
	double const * const B,
	double const * const tm,
	double const * const tn,
	double       * const tf,
	double const * const tg,
	double const * const lambda,
	double const * const Da,
	double const * const Db,
	int    const * const nobs,
	int    const * const N,
	int    const * const M,
	int    const * const ka,
	int    const * const kb,
	double const * const weight,
	int    const * const kr,
	int    const * const p,
	int const * const dl, double * dp, int * ip)
{
	double * Rw = pfdaAlloc_d(*M,&dp);
	pfda_computeResid(Rw,w,nobs, M, N, ka, B, p, tm, NULL, NULL, dl, dp);
	
	double * Rz = pfdaAlloc_d(*M,&dp);
	pfda_computeResid(Rz,z,nobs, M, N, kb, B, p, tn, NULL, NULL, dl, dp);

	double * phi = pfdaAlloc_d(*M**ka,&dp);
	dgemm_(&NoTrans,&NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);
	
	double * psi = pfdaAlloc_d(*M**kb,&dp);
	dgemm_(&NoTrans,&NoTrans, M, kb, p, &dOne, B, M, tg, p, &dzero, psi, M);
	
	double * pi = pfdaAlloc_d(*M, &dp);
	dgemv_(&NoTrans, M,p, &dOne, B, M, tm, &one, &dzero, pi, &one);
 	
	int no=0,nno=0;
	for(int i=0;i<*N;i++){
		dual_bc_w_1(w+no, ww+nno, y+no, Rw+no, Rz+no, pi+no, phi+no, psi+no, lambda, Da, Db, nobs+i, M, ka, kb, weight, kr, p, dl, dp, ip);
		no+=nobs[i];
		nno+=nobs[i]*nobs[i];
	}
}


/*! The core of the continuous/continuous model

computes the estimates and itterates until convergence

@MEMORY
	- dp = 2*p+p(ka+kb)+ka*kb+ka+kb + M + sum(nobs^2) + p^2*N + max(
	#i 	2*p^2 + M*(k+2)+ p*N + 8*p
	#1	M*(ka + kb + 3) + max(7 * max(ka,kb)^2 , 3ka +3kb + (ka + kb+ max(ka,kb))*ni)
	#2	M + M*kb + 2*kb^2
	#3	p^2+ M+ M*max(ka,kb)
	#4	M + max(ka,kb)^2 * N + p^2 + p 
	#5	ka^2 + ka*kb + 10*max(ka,kb)
	#6	4*p^2 + 2*p +  max(  ka, kb, 8 ) * p + ka*kb
	#W	M*(3+ka+kb) +   (ni*kr + kr + 2 + 5*ka + kb + 10 * k^2))
	- ip = 6*p
*/
void dual_bc_core(
	int    * const y,
	double * const B,
	double * const z,
	double * const tm,
	double * const tn,
	double * const tf,
	double * const tg,
	double * const alpha,
	double * const beta,
	double * const lambda,
	double * const Da,
	double * const Db,
	double * const sxi,
	double * const aa,
	double * const ab,
	double * const bb,
	double * const Saa,
	double * const Sab,
	double * const Sbb,
	const int * const nobs,
	const int * const N,
	const int * const M,
	const int * const ka,
	const int * const kb,
	const int * const p,
	const double * const lm,
	const double * const ln,
	const double * const lf,
	const double * const lg,
	const double * const K,
	const double * const minV,
	const int * const k0,
	const int * const kr,
	      int * const maxI,
	      double * const tol,
	const int * const dl, double * dp, int * ip)
{
	double * w = pfdaAlloc_d(*M,&dp);
	int size_ww = 0;for(int i=0;i<*N;i++)size_ww+=nobs[i]*nobs[i];
	double * ww = pfdaAlloc_d(size_ww, &dp);
	double * btb = pfdaAlloc_d(*p**p**N,&dp);
	pfda_computebtb(btb,N,B,M,p,nobs,dl);
	dual_bc_i(tm, tn, tf, tg, alpha, beta, lambda, Da, Db, sxi, aa, bb, y, z, B, nobs, btb, N, M, ka, kb, p, minV, dl, dp, ip);
	int I=0, pka = *p**ka, pkb = *p**kb, kab = *ka**kb;
	double cc=0;
	double 	sigma_xi_old=0;	
	double * tm_old     =pfdaAlloc_d(*p, &dp);
	double * tn_old     =pfdaAlloc_d(*p, &dp);
	double * tf_old     =pfdaAlloc_d(pka, &dp);
	double * tg_old     =pfdaAlloc_d(pkb, &dp);
	double * lambda_old =pfdaAlloc_d(kab, &dp);
	double * Da_old     =pfdaAlloc_d(*ka, &dp);
	double * Db_old     =pfdaAlloc_d(*kb, &dp);
	while(++I<*maxI)/*limited loop with exit*/{
		{ //setup for convergence
			sigma_xi_old = *sxi;
			dcopy_(&kab, lambda, &one, lambda_old, &one);
			dcopy_(p, tm, &one, tm_old, &one);
			dcopy_(p, tn, &one, tn_old, &one);
			dcopy_(&pka, tf, &one, tf_old, &one);
			dcopy_(&pkb, tg, &one, tg_old, &one);
			dcopy_(&kab, lambda, &one, lambda_old, &one);
			dcopy_(ka, Da, &one, Da_old, &one);
			dcopy_(kb, Db, &one, Db_old, &one);
		}
		{ // W step stuff
			double const * weight;
			if(I<*k0){
				weight = &dOne;
				dual_bc_w(w, ww, y, z, B, tm, tn, tf, tg, lambda, Da, Db, nobs, N, M, ka, kb, weight, k0, p, dl, dp, ip);
			} else {
				double cw = 10.0/(10.0+I);
				weight = &cw;
				dual_bc_w(w, ww, y, z, B, tm, tn, tf, tg, lambda, Da, Db, nobs, N, M, ka, kb, weight, kr, p, dl, dp, ip);
			}
		}
		dual_bc_1(alpha, beta, aa, ab, bb, Saa, Sab, Sbb, z, B, w, ww, tm, tn, tf, tg, lambda, Da, Db, sxi, nobs, N, M, ka, kb, p, dl, dp, ip);
		dual_bc_2(sxi, z, B, tn, tg, beta, Sbb, nobs, N, M, kb, p, dl, dp);
		dual_bc_3(tm, tn, z, B, w, tf, tg, alpha, beta, sxi, nobs, N, M, ka, kb, p, lm, ln, K, dl, dp);
		dual_bc_4(tf, tg, z, B, w, tm, tn, alpha, beta, sxi, aa, bb, nobs, N, M, ka, kb, p, lf, lg, K, btb, dl, dp, ip);
		dual_bc_5(lambda, aa, ab, N, ka, kb, dl, dp, ip);
		dual_bc_6(tf, tg, alpha, beta, lambda, Da, Db, aa, bb, N, ka, kb, p, dl, dp, ip);
		{ // Compute Convergence Criteria 
			cc = fabs(sigma_xi_old-*sxi)/MMAX(*sxi,*minV);

			daxpy_(p, &mOne, tm, &one, tm_old, &one);
			daxpy_(p, &mOne, tn, &one, tn_old, &one);
			double cctm=0;
			double cctn=0;
			for(int i=0;i<*p;i++){
				cctm+=fabs(tm_old[i])/MMAX(fabs(tm[i]),*minV);
				cctn+=fabs(tn_old[i])/MMAX(fabs(tn[i]),*minV);
			}
			cc+=cctm/(*p)+cctn/(*p);

			double cctf=0;
			{int itmp=*p**ka;daxpy_(&itmp, &mOne, tf, &one, tf_old, &one);
			for(int i=0;i<itmp;i++)cctf+=fabs(tf_old[i])/MMAX(fabs(tf[i]),*minV);
			cc+=cctf/itmp;}

			double cctg=0;
			{int itmp=*p**kb;daxpy_(&itmp, &mOne, tg, &one, tg_old, &one);
			for(int i=0;i<itmp;i++)cctg+=fabs(tg_old[i])/MMAX(fabs(tg[i]),*minV);
			cc+=cctg/itmp;}

			daxpy_(&kab, &mOne, lambda, &one, lambda_old, &one);
			double ccL=0;
			for(int i=0;i<kab;i++)ccL+=fabs(lambda_old[i])/MMAX(fabs(lambda[i]),*minV);
			cc+=ccL;

			daxpy_(ka, &mOne, Da, &one, Da_old, &one);
			double ccDa=0;
			for(int i=0;i<*ka;i++)ccDa+=cc+=fabs(Da_old[i])/MMAX(fabs(Da[i]),*minV);
			cc+=ccDa;

			daxpy_(kb, &mOne, Db, &one, Db_old, &one);
			double ccDb=0;
			for(int i=0;i<*kb;i++)ccDb+=fabs(Db_old[i])/MMAX(fabs(Db[i]),*minV);
			cc+=ccDb;

			cc/=9;
			///* Debug */
			if(checkdebug(dl,debugnum_dual_criteria_components)){
				pfda_debug_msg("Convergence criteria components:\n");
				pfda_debug_msg("sxi:  \t%5.5g\n", fabs(sigma_xi_old-*sxi));
				pfda_debug_msg("tm:   \t%5.5g\n", cctm);
				pfda_debug_msg("tn:   \t%5.5g\n", cctn);
				pfda_debug_msg("Tf:   \t%5.5g\n", cctf);
				pfda_debug_msg("Tg:   \t%5.5g\n", cctg);
				pfda_debug_msg("lambda:    \t%5.5g\n", ccL);
				pfda_debug_msg("Da:   \t%5.5g\n", ccDa);
				pfda_debug_msg("Db:   \t%5.5g\n", ccDb);
				fflush(stdout);
			}
			if(checkdebug(dl,debugnum_dual_criteria)){pfda_debug_msg("Criteria:%g\n", cc);fflush(stdout);}
		}
		if(cc < *tol){
			if(checkdebug(dl,debugnum_dual)){pfda_debug_msg("Criteria met leaving loop.%g\n");fflush(stdout);}
			break;
		}
	
	}
	*tol=cc;
	*maxI = I;	
}

