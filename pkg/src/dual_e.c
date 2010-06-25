/*! \file dual_e.c
	E step functions for the bivariate case.
*/

#include "pfda.h"


/*!  computes \f$ Sigma_eta = Db - \Lambda Da \Lambda^T \f$
@MEMORY:
	-dp length = \f$ k_a^2 + k_a *k_b \f$
*/
void pfda_gen_e_eta(
	double       * const sigma_eta,
	double const * const Db,
	double const * const Lambda,
	double const * const Da,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,
	double * dp
)
{	
	pfda_debug_cdl(debugnum_pfda_gen_e_eta){pfda_debug_msg("compute Sigma_eta:\n");fflush(stdout);}
	double * Da_matrix = pfdaAlloc_d(*ka**ka,&dp);
	addDiag(ka, Da, Da_matrix, ka);
	pfda_matrix_outer_quadratic_form(sigma_eta,Lambda,kb,kb,Da_matrix,ka, dl, dp);
	int kb2 = *kb**kb;
	dscal_(&kb2, &mOne, sigma_eta,&one); 
	addDiag(kb, Db, sigma_eta, kb);
	pfda_debug_cdl(debugnum_pfda_gen_e_eta){pfda_debug_argmat(sigma_eta,*kb,*kb);}
}

/*!  computes \f[ \Sigma_\eta = (Db - \Lambda D_\alpha \Lambda^T)^{-1} \f]
@MEMORY:
	- dp length = \f$ k_a^2 + k_a *k_b \f$ 
*/
void pfda_gen_e_eta_inv(
	double       * const sigma_eta_inv,
	double const * const Db,
	double const * const Lambda,
	double const * const Da,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,
	double * dp, int * ip
){
	pfda_gen_e_eta(sigma_eta_inv,Db,Lambda,Da,ka,kb,dl,dp);
	int sr;
	pfda_sym_inverse(sigma_eta_inv, kb, &sr, dl, dp,ip);
	if(sr){pfda_error("PFDA ERR - pfda_gen_e_eta: pfda_sym_inverse returned %i (inv_sigma_eta)",sr);}
}

/*! \f[ \zeta_{aa} = D_a^-1 + \Lambda^T Sigma_eta \Lambda + phi_i^T phi_i \f]
@MEMORY
	-dp length = \f$ k_a*k_b \f$
*/
void pfda_gen_e1_aa(
	double       * const zeta_aa,
	double const * const lambda,
	double const * const Da,
	double const * const phi,
	double const * const epsilon,
	double const * const Db,
	double const * const inv_sigma_eta,
	int    const * const M,
	int    const * const n,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,double * dp)
{ pfda_debug_step
	double inv_epsilon = epsilon?1/(*epsilon):1;
	pfda_matrix_inner_quadratic_form(zeta_aa,lambda,ka,kb,inv_sigma_eta,kb,dl,dp);
	pfda_debug_cdl(debugnum_pfda_gen_e1_aa){pfda_debug_argmat(zeta_aa,*ka,*ka);}
	dsyrk_(&Upper, &Trans, ka, n, &inv_epsilon, phi, M, &dOne, zeta_aa, ka);
	pfda_debug_cdl(debugnum_pfda_gen_e1_aa){pfda_debug_msg("zeta_aa (L^t S L+phi^t phi/epsilon):\n");pfda_debug_argmat(zeta_aa,*ka,*ka);}
	addInvDiag(ka, Da, zeta_aa, ka);
	pfda_fillsym(zeta_aa,ka,dl);
	pfda_debug_cdl(debugnum_pfda_gen_e1_aa){pfda_debug_argmat(zeta_aa,*ka,*ka);}
}

/*! \f$ \zeta_{\alpha\beta} = - \Lambda^T \Sigma_\eta^{-1} \f$
@MEMORY
	- NONE
*/
void pfda_gen_e1_ab(
	double       * const zeta_ab,
	double const * const lambda,
	double const * const inv_sigma_eta,
	int    const * const ka,
	int    const * const kb,

	int const * const dl)
{	pfda_debug_step
	dgemm_(&Trans, &NoTrans, ka, kb, kb, &mOne,  lambda, kb, inv_sigma_eta, kb, &dzero, zeta_ab, ka);	
	pfda_debug_cdl(debugnum_pfda_gen_e1_ab){pfda_debug_argmat(zeta_ab,*ka,*kb);
	                                                   pfda_debug_msg("leaving pfda_gen_e1_ab:\n");fflush(stdout);}
}

/*! \f$ \zeta_{\beta\beta} = \Sigma_\eta^{-1} +\frac{1}{\xi} \psi_i^T psi_i \f$
@MEMORY
	- NONE
*/
void pfda_gen_e1_bb(
	double       * const zeta_bb,
	double const * const inv_sigma_eta,
	double const * const psi,
	double const * const xi,
	int    const * const M,
	int    const * const n,
	int    const * const kb,
	int const * const dl)
{ pfda_debug_step
	double inv_xi = xi?1/(*xi):1;
	dsyrk_(&Upper, &Trans, kb, n, &inv_xi, psi, M, &dzero, zeta_bb, kb);
	pfda_debug_cdl(debugnum_pfda_gen_e1_bb){pfda_debug_msg("crossprod(psi):\n");pfda_debug_argmat(zeta_bb,*kb,*kb);fflush(stdout);}
	int kb2=*kb**kb;
	daxpy_(&kb2, &dOne, inv_sigma_eta, &one, zeta_bb, &one);
	pfda_fillsym(zeta_bb,kb,dl);
	pfda_debug_cdl(debugnum_pfda_gen_e1_bb){pfda_debug_argmat(zeta_bb,*kb,*kb);fflush(stdout);}
}

/*! Symetric matrices are returned in the upper portion
	\f{eqnarray*}{
		\zeta_{\alpha\alpha} &=& Da + \Lambda^T \Sigma_eta^{-1} \Lambda +\frac{1}{\sigma_\epsilon} \phi^T\phi\\
		\zeta_{\alpha\beta} &=& -\Lambda^T \Sigma_eta^{-1}\\
		\zeta_{\beta\beta} &=& \Sigma_eta^{-1} \frac{1}{\sigma_\xi} \psi^T\psi
	\f}
@MEMORY
	- dp = \f$ k_b^2 + k_a^2 + \max ( k_a*k_b  , 2*k_b^2) \f$
	- ip = \f$ k_b \f$
@callgraph
@callergraph
*/
void pfda_dual_e1(
	int    const * const M, ///< total number of observations
	int    const * const n, ///< number of subjects
	double const * const lambda, ///< \f$ \Lambda \f$
	double const * const Da, ///< \f$ D_\alpha \f$
	double const * const phi, ///< \f$ \phi=B\Theta_f \f$
	double const * const epsilon, ///< \f$ \sigma_\epsilon \f$
	int    const * const ka, ///< number of columns of \f$ phi \f$
	double const * const Db, ///< \f$ D_\beta \f$
	double const * const psi, ///< \f$ \psi = B\Theta_g \f$
	double const * const xi, ///< \f$ \sigma_\xi \f$
	int    const * const kb, ///< number of columns of \f$ \psi \f$
	double       * const zeta_aa, ///< [out] \f$ \zeta_{\alpha\alpha} \f$
	double       * const zeta_ab, ///< [out] \f$ \zeta_{\alpha\beta}  \f$
	double       * const zeta_bb, ///< [out] \f$ \zeta_{\beta\beta}  \f$
	int const * const dl, ///< debuglevel
	double * dp,  ///< double pool
	int * ip) ///< integer pool
{//code
	pfda_debug_cdl(debugnum_dual_e1){
		pfda_debug_step
		pfda_debug_arg(*M);
		pfda_debug_arg(*n);
		pfda_debug_argmat(lambda,*kb,*ka);
		pfda_debug_argmat(Da,one,*ka);
		pfda_debug_argmat2(*n,*ka,phi,M);
		pfda_debug_arg(*epsilon);
		pfda_debug_arg(*ka);
		pfda_debug_argmat(Db,one,*kb);
		pfda_debug_argmat2(*n,*kb,psi,M);
		pfda_debug_arg(*xi);
		pfda_debug_arg(*kb);
		fflush(stdout);
	} 
	int kb2=*kb**kb;
	// double inv_epsilon = epsilon?1/(*epsilon):1;   pfda_debug_cdl(debugnum_dual_e1){pfda_debug_arg(inv_epsilon);fflush(stdout);}
	// double inv_xi = xi?1/(*xi):1;             pfda_debug_cdl(debugnum_dual_e1){pfda_debug_arg(inv_xi);fflush(stdout);}

	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_msg("compute Sigma_eta:\n");fflush(stdout);}
	double * inv_sigma_eta = pfdaAlloc_d(kb2,&dp);
	pfda_gen_e_eta_inv(inv_sigma_eta,Db,lambda,Da,ka,kb,dl,dp,ip);	
	
	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_msg("compute zeta_aa:\n");fflush(stdout);}
	pfda_gen_e1_aa( zeta_aa, lambda, Da, phi, epsilon, Db, inv_sigma_eta, M, n, ka, kb, dl, dp);

	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_msg("compute zeta_ab:\n");fflush(stdout);}
	pfda_gen_e1_ab( zeta_ab, lambda, inv_sigma_eta, ka, kb, dl);
	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_argmat(zeta_ab,*ka,*kb);fflush(stdout);}
	
	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_msg("compute zeta_bb:\n");fflush(stdout);}
	pfda_gen_e1_bb( zeta_bb, inv_sigma_eta, psi, xi, M, n, kb, dl);
	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_argmat(zeta_bb,*kb,*kb);fflush(stdout);}
	
	pfda_debug_cdl(debugnum_dual_e1){pfda_debug_msg("leaving pfdaDual_e1\n");fflush(stdout);}
}

/*! Computes the components as follows:
\f{eqnarray*}{
	D&=& \begin{cases}(A-B C^{-1} B^T)^{-1} & \text{if outer=TRUE}\\
          (A-B^T C^{-1} B)^{-1} & \text{otherwise.} \end{cases}\\
	C&=&C^{-1}
\f}

debug level number = 281
@MEMORY:
	- dp	length = \f$ \max( 2 a_1^2, 2 c_1^2) \f$
	- ip	length =  \f$ \max ( a_1, c_1 ) \f$
@callgraph
@callergraph
*/
void pfda_gen_e2_1(
	int    const * const outer,///< logical determining which operation to perfom
	int    const * const a1,   ///< dimention \f$ a_1 \f$
	int    const * const c1,   ///< dimention \f$ c_1 \f$
	double const * const A,    ///< \f$ (a_1\times a_1) \f$ matrix
	double const * const B,    ///< \f$ (a_1 \times c_1) \f$ or \f$ (c_1 \times a_1) \f$ matrix depending on outer
	double       * const C,    ///<  [in/out] \f$ (c_1 \times c_1) \f$ matrix, is inverted on output.
	double       * const D,    ///<  [out] \f$(a_1 \times a_1) \f$ matrix
	int const * const dl,      ///< debuglevel
	double * dp,               ///< double pool
	int * ip)                  ///< integer pool
{/// \par Code:
	pfda_debug_step
	pfda_debug_cdl(debugnum_pfda_gen_e2_1){
		pfda_debug_msg("pfdaDual_e2_1 - \n");
		pfda_debug_msg("outer:\n%d=%s\n\n",*outer,(*outer)?("use outer"):("use inner"));
		pfda_debug_arg(*a1);
		pfda_debug_arg(*c1);
		pfda_debug_argmat(A, *a1, *a1);
		if(outer){ pfda_debug_argmat(B, *a1, *c1);} else {pfda_debug_argmat(B, *c1, *a1);}
		pfda_debug_argmat(C, *c1, *c1);
		fflush(stdout);
	}
	int sr;
	int a2=*a1**a1;

	//invert C
	double * C_inv = C;
	pfda_sym_inverse(C_inv, c1, &sr, dl, dp, ip); /// \f[ C=C^{-1} \f]
	if(sr){pfda_error("PFDA ERR - pfdaDual_e2_1: Leading minor of order %i is not positive definite (C)",sr);}
	pfda_debug_cdl(debugnum_pfda_gen_e2_1){pfda_debug_argmat(C_inv, *c1, *c1); fflush(stdout);}
	
	/*! \f[ 
		D = \begin{cases} B C^{-1} B^T	& \text{if outer=TRUE} \\ B^T C^{-1} B	& \text{otherwise.} \end{cases}
	\f] */
	if(*outer){ pfda_matrix_outer_quadratic_form(D, B, a1, a1, C_inv, c1, dl, dp); 
	} else {    pfda_matrix_inner_quadratic_form(D, B, a1, c1, C_inv, c1, dl, dp); }
	pfda_debug_cdl(debugnum_pfda_gen_e2_1){
		if(*outer){pfda_debug_msg("B(C_inv)B^t:"); pfda_debug_argmat(D, *a1, *a1); 
		}else{     pfda_debug_msg("B^t(C_inv)B:"); pfda_debug_argmat(C_inv, *a1, *a1);}
		fflush(stdout);
	}
	dscal_(&a2, &mOne, D,&one); /// \f[ D=-D \f]
	daxpy_(&a2, &dOne, A, &one, D, &one);  /// \f[ D = D + A \f]
	pfda_debug_cdl(debugnum_pfda_gen_e2_1){pfda_debug_msg("A-i/o_q_form(B,C): "); pfda_debug_argmat(D, *a1, *a1); fflush(stdout);}

	pfda_sym_inverse(D, a1, &sr, dl, dp, ip);  /// \f[ D=D^{-1} \f]
	if(sr){pfda_error("PFDA ERR - pfdaDual_e2_1: Leading minor of order %i is not positive definite (D)",sr);}
	pfda_debug_cdl(debugnum_pfda_gen_e2_1){pfda_debug_msg("D=solve(A-i/o_q_form(B,C)): \n"); pfda_debug_argmat(D, *a1, *a1); fflush(stdout);}
	pfda_debug_cdl(debugnum_pfda_gen_e2_1){pfda_debug_msg("leaving pfda_gen_e2_1\n"); fflush(stdout);}
}

/*! \f$(\zeta_{\alpha\alpha} - \zeta_{\alpha\beta} \zeta_{\beta\beta}^{-1} \zeta_{\alpha\beta}^T)^{-1} \f$
@MEMORY
 	- dp	length = \f$ \max( 2 k_a^2, 2 k_b^2) \f$
	- ip	length =  \f$ \max ( k_a, k_b ) \f$
*/
void pfda_gen_e2_aa(
	double       * const sigma_aa,
	double const * const zeta_aa,
	double const * const zeta_ab,
	double       * const zeta_bb,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,

	double * dp, int * ip)
{
	pfda_gen_e2_1(&one,  ka, kb, zeta_aa, zeta_ab, zeta_bb, sigma_aa, dl, dp, ip);
}

/*! \f$ (\zeta_{\alpha\alpha} - \zeta_{\alpha\beta} (\zeta_{\beta\beta})^{-1} \zeta_{\alpha\beta}^T)^{-1} \f$
@MEMORY
 	- dp	length = \f$ \max( 2 k_a^2, 2 k_b^2) \f$ 
	- ip	length =  \f$ \max( k_a, k_b ) \f$ 
*/
void pfda_gen_e2_bb(
	double * const sigma_bb,
	double       * const zeta_aa,
	double const * const zeta_ab,
	double       * const zeta_bb,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,

	double * dp, int * ip)
{
	pfda_gen_e2_1(&zero, kb, ka, zeta_bb, zeta_ab, zeta_aa, sigma_bb, dl, dp, ip);
}

/*! \f$ \Sigma _{\alpha\beta} = - \Sigma_{\alpha\alpha} \zeta_{\alpha\beta} (\zeta_{\beta\beta})^{-1} \f$
@MEMORY
	- dp length  = \f$ k_a*k_b \f$
*/
void pfda_gen_e2_ab(
	double       * const sigma_ab,
	double const * const sigma_aa,
	double const * const zeta_ab,
	double const * const zeta_bb_inv,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,
	double * dp)
{
	pfda_debug_cdl(debugnum_pfda_gen_e2_ab){pfda_debug_msg("entering pfda_gen_e2_ab\n"); fflush(stdout);}
	double * tmp = pfdaAlloc_d(*ka**kb, &dp);
	dgemm_(&NoTrans, &NoTrans, ka,kb,ka, &dOne,  sigma_aa, ka, zeta_ab, ka, &dzero, tmp, ka);
	pfda_debug_cdl(debugnum_pfda_gen_e2_ab){pfda_debug_msg("sigma_aa zeta_ab:\n"); pfda_debug_argmat(tmp, *ka, *kb); fflush(stdout);}
	dsymm_(&Right, &Upper, ka, kb, &mOne, zeta_bb_inv, kb, tmp, ka, &dzero, sigma_ab, ka);	
	pfda_debug_cdl(debugnum_pfda_gen_e2_ab){pfda_debug_msg("leaving pfda_gen_e2_ab\n"); fflush(stdout);}
}

/*! Computes Variance matrices
	\f{eqnarray*}{
		\Sigma_{\alpha\alpha} &=& (\zeta_aa-\zeta_ab \zeta_bb^{-1} \zeta_ab^T)^{-1}\\
		\Sigma_{\alpha\beta}  &=& -\sigma_aa\zeta_ab \zeta_bb^{-1}\\
		\Sigma_{\beta\beta}   &=& (\zeta_bb-\zeta_ab^T \zeta_aa^{-1} \zeta_ab)^{-1}
	\f}
@MEMORY:
	- dp	length = \f$ k_b^2 + \max ( 2k_a^2, 2*k_b^2 ) \f$   use \f$ 3* \max(k_a,k_b)^2 \f$ 
	- ip	length = \f$ \max ( k_a, k_b) \f$ 
*/
void pfdaDual_e2(
	int    const * const ka,      ///< \f$ k_\alpha \f$ 
	int    const * const kb,      ///< \f$ k_\beta \f$
	double       * const zeta_aa, ///< \f$ \zeta_{\alpha\alpha} \f$
	double const * const zeta_ab, ///< \f$ \zeta_{\alpha\beta} \f$
	double       * const zeta_bb,///< \f$ \zeta_{\beta\beta} \f$
	double       * const sigma_aa,///< \f$ \Sigma_{\alpha\alpha} \f$
	double       * const sigma_ab,///< \f$ \Sigma_{\alpha\beta} \f$
	double       * const sigma_bb,///< \f$ \Sigma_{\beta\beta}  \f$
	int const * const dl, ///< debug level
	double * dp, ///< double pool
	int * ip	) ///< integer pool
{// \par Code:
	pfda_debug_step;
	pfda_debug_cdl(debugnum_dual_e2){
		pfda_debug_msg("pfdaDual_e2 - \n");
		pfda_debug_arg(*ka);
		pfda_debug_arg(*kb);
		pfda_debug_argmat(zeta_aa, *ka, *ka);
		pfda_debug_argmat(zeta_ab, *ka, *kb);
		pfda_debug_argmat(zeta_bb, *kb, *kb);
		fflush(stdout);
	}
 	int kb2=*kb**kb;
	double * zeta_bb_inv = pfdaAlloc_d(kb2,&dp);
	
	dcopy_(&kb2,zeta_bb, &one, zeta_bb_inv, &one);
	pfda_debug_cdl(debugnum_dual_e2){pfda_debug_argmat(zeta_bb_inv, *kb, *kb); fflush(stdout);}
	pfda_gen_e2_aa( sigma_aa, zeta_aa, zeta_ab, zeta_bb_inv, ka, kb, dl, dp, ip);
	pfda_debug_cdl(debugnum_dual_e2){pfda_debug_argmat(zeta_bb_inv, *kb, *kb); fflush(stdout);}
	pfda_gen_e2_bb( sigma_bb, zeta_aa, zeta_ab, zeta_bb, ka, kb, dl, dp, ip);
	pfda_debug_cdl(debugnum_dual_e2){pfda_debug_argmat(sigma_aa, *ka, *ka); fflush(stdout);}
	pfda_debug_cdl(debugnum_dual_e2){pfda_debug_argmat(sigma_bb, *kb, *kb); fflush(stdout);}
	
	pfda_gen_e2_ab( sigma_ab, sigma_aa, zeta_ab, zeta_bb_inv, ka, kb, dl, dp);
	pfda_debug_cdl(debugnum_dual_e2){
		pfda_debug_argmat(sigma_ab, *ka, *kb); 
		pfda_debug_msg("leaving pfdaDual_e2\n");fflush(stdout);}
}

/*! \f[ \mu = \frac{1}{\sigma_\epsilon} * \Sigma_{AA} \psi u + \frac{1}{\sigma_\xi} * \Sigma_{AB} \phi v \f]
 *	@MEMORY
 *		-dp	\f$ 2 a*n \f$
 */
void pfda_gen_e3_1(
	double       * const mu,     ///< [out] \f$ \mu \f$
	int    const * const ldmu,   ///< increment of \f$ \mu \f$
	double const * const epsilon,///< \f$ \sigma_\epsilon \f$
	double const * const AA,     ///< \f$ \Sigma_{AA} \f$ matrix of dimention \f$(a \times a)\f$
	double const * const phi,    ///< \f$ \phi \f$
	int    const * const ldphi,  ///< Total number of rows of \f$ \phi \f$
	double const * const u,      ///< vector of length nu
	int    const * const nu,     ///< number of rows of u
	int    const * const a,      ///< dimention of AA
	double const * const xi,     ///< \f$ \sigma_\xi \f$
	double const * const AB,     ///< \f$ \Sigma_{AB} \f$ matrix of dimention \f$ (a \times b) \f$.
	int    const * const transAB,///< logical indicating if AB should be transposed of not.
	double const * const psi,    ///< \f$ \psi \f$
	int    const * const ldpsi,  ///< Total number of rows of \f$ \psi \f$
	double const * const v,      ///< vector of length nv
	int    const * const nv,     ///< number of observations for v
	int    const * const b,      ///< number of columns of AB.
	int const * const dl,        ///< debug level
	double * dp)                 ///< double pool
{//code
	pfda_debug_cdl(debugnum_pfda_gen_e3_1){
		pfda_debug_step
		pfda_debug_arg(*epsilon);
		pfda_debug_argmat(AA, *a, *a);
		pfda_debug_argmat2(*nu,*a, phi,ldphi);
		pfda_debug_arg(*ldphi);
		pfda_debug_argmat(u, one, *nu);
		pfda_debug_arg(*nu);
		pfda_debug_arg(*a);
		pfda_debug_argmat(xi, one, one);
		pfda_debug_argmat(AB, *a, *b);
		pfda_debug_argmat2(*nv,*b, psi,ldpsi);
		pfda_debug_arg(*ldpsi);
		pfda_debug_argmat(v, one, *nv);
		pfda_debug_arg(*nv);
		pfda_debug_arg(*b);
		fflush(stdout);
	}
	double * Xa = pfdaAlloc_d(*a**nu,&dp);
	double * Xb = pfdaAlloc_d(*a**nv,&dp);
	double inv_epsilon = epsilon ? 1/(*epsilon) : 1;   pfda_debug_cdl(debugnum_pfda_gen_e3_1){pfda_debug_arg(inv_epsilon);fflush(stdout);}
	double inv_xi      = xi      ? 1/(*xi)      : 1;             pfda_debug_cdl(debugnum_pfda_gen_e3_1){pfda_debug_arg(inv_xi);fflush(stdout);}
	
	dgemm_(&NoTrans, &Trans, a, nu, a, &inv_epsilon, AA, a, phi, ldphi, &dzero, Xa, a);
	dgemm_((*transAB)?(&Trans):(&NoTrans), &Trans, a, nv, b, &inv_xi, AB, (*transAB)?b:a, psi, ldpsi, &dzero, Xb, a);
	pfda_debug_cdl(debugnum_pfda_gen_e3_1){pfda_debug_argmat(Xa, *a, *nu); fflush(stdout);}
	pfda_debug_cdl(debugnum_pfda_gen_e3_1){pfda_debug_argmat(Xb, *a, *nv); fflush(stdout);}
	
	dgemv_(&NoTrans, a, nu, &dOne, Xa, a, u, &one, &dzero, mu, ldmu); 
	pfda_debug_cdl(debugnum_pfda_gen_e3_1){pfda_debug_argmat2(one, *a,mu,ldmu); fflush(stdout);}
	dgemv_(&NoTrans, a, nv, &dOne, Xb, a, v, &one, &dOne , mu, ldmu);
	pfda_debug_cdl(debugnum_pfda_gen_e3_1){
		pfda_debug_argmat2(one, *a,mu,ldmu);
		pfda_debug_msg("leaving pfda_gen_e3_1\n");fflush(stdout);}
}

/*! \f[ \mu  = \frac{1}{\sigma_\epsilon} \Sigma_{AA} \psi u + \frac{1}{\sigma_\xi} \Sigma_{AB} \phi v \f]
 *	@MEMORY
 *		- dp	\f$ 2 a*n \f$
 */
void pfda_dual_e3_1(
	double       * const mu,     ///< [out] \f$ \mu \f$
	int    const * const ldmu,   ///< incrament for \f$ \mu \f$
	int    const * const M,      ///< total number of observations
	int    const * const n,      ///< number of observations for subject
	double const * const epsilon,///< \f$ \sigma_\epsilon \f$
	double const * const AA,     ///< \f$ \Sigma_{AA} \f$ symetric matrix \f$(a \times a)\f$
	double const * const phi,    ///< \f$ \phi \f$
	double const * const u,      ///< vector of length a 
	int    const * const a,      ///< length of u
	double const * const xi,     ///< \f$ \sigma_\xi \f$
	double const * const AB,     ///< \f$ \Sigma_{AB} \f$ matrix \f$ (a \times b) \f$
	int    const * const transAB,///< logical - does AB need to be transposed.
	double const * const psi,    ///< \f$ \psi \f$
	double const * const v,      ///< vector of length b
	int    const * const b,      ///< number of rows of b
	int const * const dl,        ///< debug level
	double * dp)                 ///< double pool
{//code
	pfda_gen_e3_1( mu, ldmu, epsilon, AA, phi, M, u, n, a, xi, AB, transAB, psi, M, v,  n, b, dl, dp);
}

/*! Computes the variance matrices of alpha and beta given the data.
	\f{eqnarray*}{
		\zeta_{\alpha\alpha} &=& Da + \Lambda^T \Sigma_eta^{-1} \Lambda +\frac{1}{\sigma_\epsilon} \phi^T\phi\\
		\zeta_{\alpha\beta} &=& -\Lambda^T \Sigma_eta^{-1}\\
		\zeta_{\beta\beta} &=& \Sigma_eta^{-1} \frac{1}{\sigma_\xi} \psi^T\psi \\
		\Sigma_{\alpha\alpha} &=& (\zeta_{\alpha\alpha}-\zeta_{\alpha\beta} \zeta_{\beta\beta}^{-1} \zeta_{\alpha\beta}^T)^{-1}\\
		\Sigma_{\alpha\beta}  &=& -\Sigma_{\alpha\alpha}\zeta_{\alpha\beta} \zeta_{\beta\beta}^{-1}\\
		\Sigma_{\beta\beta}   &=& (\zeta_{\beta\beta}-\zeta_{\alpha\beta}^T \zeta_{\alpha\alpha}^{-1} \zeta_{\alpha\beta})^{-1}
	\f}
@MEMORY
	- dp \f$ k_a^2 + k_b^2 + k_a*k_b + \max( \f$
		-# \f$ k_b^2 + k_a^2 + \max ( k_a*k_b  , 2 k_b^2)  \f$ 
		-#  \f$ k_b^2 + \max ( 2 k_a^2, 2 k_b^2 )  \f$   use  \f$ 3 \max(k_a,k_b)^2) \f$ 
		-  \f$ 7 \max(k_a, k_b)^2 \f$  is ample to cover it.
	- ip \f$ \max(k_a, k_b) \f$
@callgraph
@callergraph
*/
void dual_gen_sigmas(
				double * const Saa,   ///< \f$ \Sigma_{\alpha\alpha} \f$
				double * const Sab,   ///< \f$ \Sigma_{\alpha\beta} \f$
				double * const Sbb,   ///< \f$ \Sigma_{\beta\beta} \f$
	const double * const phi,   ///< \f$ \phi = B\Theta_f \f$
	const double * const psi,   ///< \f$ \psi=B\Theta_g \f$
	const double * const lambda,///< \f$  \Lambda \f$
	const double * const Da,    ///< \f$ D_\alpha \f$
	const double * const Db,    ///< \f$ D_\beta \f$
	const double * const sep,   ///< \f$ \sigma_\epsilon \f$
	const double * const sxi,   ///< \f$ \sigma_\xi \f$
	const int    * const M,     ///< Total number of observations and the nubmer of rows of \f$\phi\f$ and \f$\psi\f$.
	const int    * const ni,    ///< \f$ n_i \f$ the number of observations for the current subject
 	const int    * const ka,    ///< number of principal components for Y
	const int    * const kb,    ///< number of principal components for Z
	const int * const dl,       ///< debug level
	double*dp,                  ///< double pool
	int * ip)                   ///< integer pool
{ pfda_debug_step;
	double * zeta_aa = pfdaAlloc_d(*ka**ka,&dp);
	double * zeta_ab = pfdaAlloc_d(*ka**kb,&dp);
	double * zeta_bb = pfdaAlloc_d(*kb**kb,&dp);
	pfda_dual_e1(M, ni, lambda, Da,  phi,   sep, ka, Db, psi , sxi, kb, zeta_aa, zeta_ab, zeta_bb, dl,dp, ip);
	pfdaDual_e2(ka,kb,zeta_aa,zeta_ab,zeta_bb,Saa,Sab,Sbb,dl,dp,ip);
	pfda_debug_lstep;
}

/*! Update \f$ \alpha \f$ and \f$ \beta \f$
	\f{eqnarray*}{
		\alpha = \frac{1}{\sigma_\epsilon} \Sigma_{ \alpha\alpha } \phi u + \frac{1} {\sigma_\xi} \Sigma_{\alpha\beta } \phi v   \\
		\beta  = \frac{1}{\sigma_\epsilon} * \Sigma_{ \alpha\beta }^T \phi u + \frac{1} {\sigma_\xi} \Sigma_{\beta\beta} \psi v
	\f}
	@MEMORY
		-dp	\f$ 2*n * \max(a,b) \f$
*/
void pfda_dual_e3(
	double       * const alpha,   ///<  \f$ \alpha \f$
	double       * const beta,    ///<  \f$ \beta \f$
	double const * const sigma_aa,///<  \f$ \Sigma_{\alpha\alpha} \f$
	double const * const sigma_ab,///<  \f$ \Sigma_{\alpha\beta} \f$
	double const * const sigma_bb,///<  \f$ \Sigma_{\beta\beta} \f$
	double const * const phi,     ///<  \f$ phi \f$
	double const * const psi,     ///<  \f$ \psi \f$
	double const * const epsilon, ///<  \f$ \sigma_\epsilon \f$
	double const * const xi,      ///<  \f$ \sigma_\xi \f$
	double const * const u,       ///<  \f$ u \f$
	double const * const v,       ///<  \f$ v \f$
	int    const * const ka,      ///<  \f$ k_a \f$
	int    const * const kb,      ///<  \f$ k_b \f$
	int    const * const n,       ///<  \f$ n_i \f$ numer of obervations for this subject ( nrow or psi and phi)
	int    const * const N,       ///<  Number of subjects
	int    const * const M,       ///<  Total number of observations
	int const * const dl,         ///<  debug level
	double * dp)                  ///<  double pool
{//code
	pfda_debug_cdl(debugnum_dual_e3){pfda_debug_msg("Entering pfdaDual_e3\n");fflush(stdout);}
	pfda_dual_e3_1(alpha,N,M,n,epsilon,sigma_aa,phi,u,ka,xi     ,sigma_ab,&zero,psi,v,kb,dl,dp);
	pfda_dual_e3_1(beta ,N,M,n,xi     ,sigma_bb,psi,v,kb,epsilon,sigma_ab,&one ,phi,u,ka,dl,dp);
	pfda_debug_cdl(debugnum_dual_e3){
		pfda_debug_argmat2(one, *ka,alpha, N);
		pfda_debug_argmat2(one, *kb,beta, N);
		pfda_debug_msg("leaving pfdaDual_e3 \n");
		fflush(stdout);
	}
}

/*! computes the Estep of the Paired fucntional data analysis problem to Zhou, et al.
@MEMORY
	-dp	\f$ 2 M + M k_a + M k_b+k_a^2 + k_a k_b + k_b^2 + \max( \f$
		-#	step 1.	\f$ k_b^2 + k_a^2 + \max ( k_a*k_b  , 2*k_b^2)  < 4\max(k_a,k_b)^2 \f$
		-#	step 2.	\f$ k_b^2 + max ( 2 k_a^2, 2 k_b^2 )  < 3\max(k_a,k_b)^2 \f$
		-#	step 3.	\f$ 2 *\max_i(n_i) * max(a,b) \f$
			= \f$ \max( 2*\max_i(n_i)*\max(a,b), 4*\max(k_a,k_b)^2) \f$  might give a little extra, the but the extra should be really small provided ka and kb are small and/or close to each other.
	-ip		\f$ \max ( k_a, k_b) \f$
*/
void pfdaDual_e(
	/* input Values */
	double const * const y,   ///< [in] \f$ y \f$ response Y
	double const * const z,   ///< [in]  \f$ z \f$ response Z
	int    const * const nobs,///< [in]  \f$ \{n_1,\dots,n_N\} \f$ vector of the number of observations for each 
	int    const * const M,   ///< [in]  Total number of observations
	int    const * const N,   ///< [in]  Number of subjects
	int    const * const ka,  ///< [in]  \f$ k_a \f$ the number of principal components for Y.
	int    const * const kb,  ///< [in]  \f$ k_b \f$ the number of principal components for Z.
	const double * const B,   ///< [in]  \f$ B \f$ The matrix of evaluated basis functionsat the times points.
	int    const * const p,   ///< [in]  Number of columns of B
	/* State Values */
	double const * const tm,      ///< [in]  \f$ \theta_\mu \f$
	double const * const tn,      ///< [in]  \f$ \theta_\nu \f$
	double const * const tf,      ///< [in]  \f$ \Theta_f \f$
	double const * const tg,      ///< [in]  \f$ \Theta_g \f$
	double const * const Da,      ///< [in]  \f$ D_\alpha \f$
	double const * const Db,      ///< [in]  \f$ D_\beta \f$
	double const * const Lambda,  ///< [in]  \f$ \Lambda \f$
	double const * const seps,    ///< [in]  \f$ \sigma_\epsilon \f$
	double const * const sxi,     ///< [in]  \f$ \sigma_\xi \f$
	double       * const alpha,   ///< [out]  matrix of \f$ \alpha_i \f$ values.
	double       * const Beta,    ///< [out]  matrix of \f$ \beta_i \f$ values.
 	double       * const Sigma_aa,///< [out]  \f$ \Sigma_{\alpha\alpha} \f$
	double       * const Sigma_ab,///< [out]  \f$ \Sigma_{\alpha\beta} \f$
	double       * const Sigma_bb,///< [out]  \f$ \Sigma_{\beta\beta} \f$
	/* Control Values */
	int const * const dl,///< debug level
	double * dp,         ///< double pool
	int * ip)            ///< integer pool
{ /* E Step */
 	{ //debug
	pfda_debug_cdl(debugnum_dual_e){pfda_debug_msg("Entering pfdaDual_e\n");fflush(stdout);}
	pfda_debug_cdl(debugnum_dual_e){
			pfda_debug_msg("INPUTS\n----------------------------------------\n");
			pfda_debug_argmat(y, 1, *M);
			pfda_debug_argmat(z, 1, *M);
			pfda_debug_arg(*M);
			pfda_debug_arg(*N );
			pfda_debug_arg(*ka);
			pfda_debug_arg(*kb);
			pfda_debug_arg(*p);
			pfda_debug_argmati(nobs, one, *N);
			pfda_debug_argvec(tm,p);
			pfda_debug_argvec(tn,p);
			pfda_debug_argmat(tf, *p, *ka);
			pfda_debug_argmat(tg, *p, *kb);
			pfda_debug_argvec(Da,ka);
			pfda_debug_argvec(Db,kb);
			pfda_debug_argmat(Lambda, *kb, *ka);
			pfda_debug_arg(*seps); 
			pfda_debug_arg(*sxi); 
			pfda_debug_argmat(alpha, *N, *ka);
			pfda_debug_argmat(Beta, *N, *kb);
		}
	pfda_debug_cdl(debugnum_dual_e_inputs_large){
			pfda_debug_argmat(Sigma_aa, *N, *ka**ka);
			pfda_debug_argmat(Sigma_aa, *N, *ka**ka);
			pfda_debug_argmat(Sigma_ab, *N, *ka**kb);
			pfda_debug_argmat(B, *M, *p);
			pfda_debug_msg("END INPUTS\n----------------------------------------\n");
		}
	}  
	{ //code
	int obsnum;
	double * Ry=pfdaAlloc_d(*M,&dp);
	double * Rz=pfdaAlloc_d(*M,&dp);
	double * Ryi=Ry;
	double * Rzi=Rz;
	double * phi=pfdaAlloc_d(*M**ka,&dp);
	double * psi=pfdaAlloc_d(*M**kb,&dp);
	double * phii=phi;
	double * psii=psi;
	double * zeta_aa=pfdaAlloc_d(*ka**ka,&dp);
	double * zeta_ab=pfdaAlloc_d(*ka**kb,&dp);
	double * zeta_bb=pfdaAlloc_d(*kb**kb,&dp);
	
	if(y)dcopy_(M, y, &one, Ry, &one);
	if(z)dcopy_(M, z, &one, Rz, &one);

	dgemv_(&NoTrans, M, p, &mOne, B, M, tm, &one, &dOne, Ry, &one); //Ry -= B*tm  (y-B*tm)
	dgemv_(&NoTrans, M, p, &mOne, B, M, tn, &one, &dOne, Rz, &one); //Rz -= B*tn  (z-B*tn)
	
	pfda_debug_cdl(debugnum_dual_e){pfda_debug_argvec(Ry,M); fflush(stdout);}
	pfda_debug_cdl(debugnum_dual_e){pfda_debug_argvec(Rz,M); fflush(stdout);}

	dgemm_(&NoTrans, &NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);//  phi = B*tf
	dgemm_(&NoTrans, &NoTrans, M, kb, p, &dOne, B, M, tg, p, &dzero, psi, M);//  psi = B*tg
	
	pfda_debug_cdl(debugnum_dual_e){pfda_debug_argmat(phi, *M, *ka); fflush(stdout);}
	pfda_debug_cdl(debugnum_dual_e){pfda_debug_argmat(psi, *M, *kb); fflush(stdout);}
		
	for(obsnum=0;obsnum<*N;obsnum++){
		pfda_debug_cdl(debugnum_dual_e){pfda_debug_msg("observation: %d ... ",obsnum+1);fflush(stdout);}
		pfda_dual_e1(M,	nobs+obsnum, Lambda, 
				Da, phii, seps, ka,	
				Db, psii, sxi , kb, 
				zeta_aa, zeta_ab, zeta_bb,
				dl,dp,ip);
		pfdaDual_e2(ka,kb,
				zeta_aa,zeta_ab,zeta_bb,
				Sigma_aa+obsnum**ka**ka,Sigma_ab+obsnum**ka**kb,Sigma_bb+obsnum**kb**kb,
				dl,dp,ip);
		pfda_dual_e3(alpha+obsnum,Beta+obsnum,
				Sigma_aa+obsnum**ka**ka,Sigma_ab+obsnum**ka**kb,Sigma_bb+obsnum**kb**kb,
				phii,psii,
				seps,sxi,
				Ryi,Rzi,
				ka,kb,
				nobs+obsnum,N,M,
				dl,dp);
		psii+=nobs[obsnum];
		phii+=nobs[obsnum];
		Ryi +=nobs[obsnum];
		Rzi +=nobs[obsnum];
	}
	pfda_debug_cdl(debugnum_dual_e){pfda_debug_msg("Leaving pfdaDual_e\n");fflush(stdout);}
	}
}


