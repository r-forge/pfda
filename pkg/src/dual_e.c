/*! \file dual_e.c
	E step functions for the bivariate case.
*/

#include "pfda.h"


/*!  computes Sigma_eta = Db - \Lambda Da \Lambda\Trans
MEMORY:
	dp length = ka^2 + ka *kb
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
	if(checkdebug(dl,debugnum_pfda_gen_e_eta)){pfda_debug_msg("compute Sigma_eta:\n");fflush(stdout);}
	double * Da_matrix = pfdaAlloc_d(*ka**ka,&dp);
	addDiag(ka, Da, Da_matrix, ka);
	pfda_matrix_outer_quadratic_form(sigma_eta,Lambda,kb,kb,Da_matrix,ka, dl, dp);
	pfdaFree(Da_matrix);
	int kb2 = *kb**kb;
	dscal_(&kb2, &mOne, sigma_eta,&one); 
	addDiag(kb, Db, sigma_eta, kb);
	if(checkdebug(dl,debugnum_pfda_gen_e_eta)){pfda_debug_msg("Sigma_eta:\n");printmat(sigma_eta,*kb,*kb);fflush(stdout);}
}

/*!  computes \f[Sigma_eta = (Db - \Lambda Da \Lambda\Trans)\inv\f]
MEMORY:
	dp length = ka^2 + ka *kb
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

/* \f[\zeta_{aa} = D_a^-1 + \Lambda\trans Sigma_eta \Lambda + phii\trans phii\f]
MEMORY:
	dp length = ka*kb
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
{
	if(checkdebug(dl,debugnum_pfda_gen_e1_aa)){pfda_debug_msg("entering pfda_gen_e1_aa\n");fflush(stdout);}
	double inv_epsilon = epsilon?1/(*epsilon):1;
	pfda_matrix_inner_quadratic_form(zeta_aa,lambda,ka,kb,inv_sigma_eta,kb,dl,dp);
	if(checkdebug(dl,debugnum_pfda_gen_e1_aa)){pfda_debug_msg("zeta_aa (L^tSL):\n");printmat(zeta_aa,*ka,*ka);fflush(stdout);}
	dsyrk_(&Upper, &Trans, ka, n, &inv_epsilon, phi, M, &dOne, zeta_aa, ka);
	if(checkdebug(dl,debugnum_pfda_gen_e1_aa)){pfda_debug_msg("zeta_aa (L^t S L+phi^t phi/epsilon):\n");printmat(zeta_aa,*ka,*ka);fflush(stdout);}
	addInvDiag(ka, Da, zeta_aa, ka);
	pfda_fillsym(zeta_aa,ka,dl);
	if(checkdebug(dl,debugnum_pfda_gen_e1_aa)){pfda_debug_msg("zeta_aa:\n");printmat(zeta_aa,*ka,*ka);fflush(stdout);}
}

void pfda_gen_e1_ab(
	double       * const zeta_ab,
	double const * const lambda,
	double const * const inv_sigma_eta,
	int    const * const ka,
	int    const * const kb,

	int const * const dl)
/*DESCRIPTION
	zeta_ab = - Lambda\trans(Sigma_eta)^-1
MEMORY
	NONE
*/
{	
	if(checkdebug(dl,debugnum_pfda_gen_e1_ab)){pfda_debug_msg("entering pfda_gen_e1_ab:\n");fflush(stdout);}
	dgemm_(&Trans, &NoTrans, ka, kb, kb, &mOne,  lambda, kb, inv_sigma_eta, kb, &dzero, zeta_ab, ka);	
	if(checkdebug(dl,debugnum_pfda_gen_e1_ab)){pfda_debug_msg("zetaab:\n");printmat(zeta_ab,*ka,*kb);
	                                                   pfda_debug_msg("leaving pfda_gen_e1_ab:\n");fflush(stdout);}
}

void pfda_gen_e1_bb(
	double       * const zeta_bb,
	double const * const inv_sigma_eta,
	double const * const psi,
	double const * const xi,
	int    const * const M,
	int    const * const n,
	int    const * const kb,
	int const * const dl)
/* DESCRIPTION
	zeta_bb = sigma_eta^-1 +1/xi * psii\trans psii
MEMORY
	NONE
*/
{
	double inv_xi = xi?1/(*xi):1;
	if(checkdebug(dl,debugnum_pfda_gen_e1_bb)){pfda_debug_msg("entering pfda_gen_e1_bb:\n");fflush(stdout);}
	dsyrk_(&Upper, &Trans, kb, n, &inv_xi, psi, M, &dzero, zeta_bb, kb);
	if(checkdebug(dl,debugnum_pfda_gen_e1_bb)){pfda_debug_msg("crossprod(psi):\n");printmat(zeta_bb,*kb,*kb);fflush(stdout);}
	int kb2=*kb**kb;
	daxpy_(&kb2, &dOne, inv_sigma_eta, &one, zeta_bb, &one);
	pfda_fillsym(zeta_bb,kb,dl);
	if(checkdebug(dl,debugnum_pfda_gen_e1_bb)){pfda_debug_msg("zeta_bb:\n");printmat(zeta_bb,*kb,*kb);fflush(stdout);}
}

void pfda_dual_e1(
	int    const * const M,
	int    const * const n,
	double const * const lambda,
	double const * const Da,
	double const * const phi,
	double const * const epsilon,
	int    const * const ka,
	double const * const Db,
	double const * const psi,
	double const * const xi,
	int    const * const kb,
	double       * const zeta_aa,
	double       * const zeta_ab,
	double       * const zeta_bb,
	int const * const dl,double * dp, int * ip)
/* DESCRIPTION
INPUTS:
	M - for this purpose the number of rows for psi and phi
	n - the number of observations for this subject.  The number of rows in psi and phi
	lambda (kb x ka)  - represent the covariance of alpha and beta
	Da (ka)  - represents vaiance of alpha
	phi (n x ka) - represents Theta_f %*% B_i
	epsilon - sigma_epsilon
	ka - number of principle components for y
	Db (kb)  - represents vaiance of beta
	psi (n x kb) - represents Theta_g %*% B_i
	xi - sigma_xi
	kb - number of principle components for z
OUTPUT:
	symetric matrices are returned in the upper portion
	\[ %tex description
		\zeta_aa <- Da + \Lambda\trans \Sigma_eta^{-1} \Lambda +\frac{1}{\sigma_\epsilon} \phi\trans\phi
		\zeta_ab <- -\Lambda\trans \Sigma_eta^{-1}
		\zeta_bb <- \Sigma_eta^{\1} \frac{1}{\sigma_\xi} \psi\trans\psi
	\]
MEMORY:
	dp - kb2 + ka2 + max ( ka*kb  , 2*kb2)
	ip - kb
*/
{//code
	if(checkdebug(dl,debugnum_dual_e1)){
		pfda_debug_msg("pfdaDual_e1 - \n");
		pfda_debug_msg("M:%d\n\n",*M);
		pfda_debug_msg("n:%d\n\n",*n);
		pfda_debug_msg("lambda:\n");printmat(lambda,*kb,*ka);
		pfda_debug_msg("Da:\n");printmat(Da,one,*ka);
		pfda_debug_msg("phi:\n");printmat2(*n,*ka,phi,M);
		pfda_debug_msg("epsilon:\n");printmat(epsilon,one,one);
		pfda_debug_msg("ka:%d\n\n",*ka);
		pfda_debug_msg("Db:\n");printmat(Db,one,*kb);
		pfda_debug_msg("psi:\n");printmat2(*n,*kb,psi,M);
		pfda_debug_msg("xi:\n");printmat(xi,one,one);
		pfda_debug_msg("kb:%d\n\n",*kb);
		fflush(stdout);
	} 
	int kb2=*kb**kb;
	double inv_epsilon = epsilon?1/(*epsilon):1;   if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("inv_epsilon:\n%9.8g\n\n",inv_epsilon);fflush(stdout);}
	double inv_xi = xi?1/(*xi):1;             if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("inv_xi:\n%9.8g\n\n",inv_xi);fflush(stdout);}

	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("compute Sigma_eta:\n");fflush(stdout);}
	double * inv_sigma_eta = pfdaAlloc_d(kb2,&dp);
	pfda_gen_e_eta_inv(inv_sigma_eta,Db,lambda,Da,ka,kb,dl,dp,ip);	
	
	// if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("compute inv_sigma_eta:\n");fflush(stdout);}
	// int sr;
	// pfda_sym_inverse(sigma_eta, kb, &sr, dl, dp,ip);
	// double * inv_sigma_eta = sigma_eta; //Aliasing // pfdaFree(sigma_eta);
	// if(sr){pfda_error("PFDA ERR - pfdaDual_e1_1: pfda_sym_inverse returned %i (inv_sigma_eta)",sr);}
	// if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("inv_sigma_eta:\n");printmat(inv_sigma_eta,*kb,*kb);fflush(stdout);}
	
	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("compute zeta_aa:\n");fflush(stdout);}
	pfda_gen_e1_aa( zeta_aa, lambda, Da, phi, epsilon, Db, inv_sigma_eta, M, n, ka, kb, dl, dp);

	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("compute zeta_ab:\n");fflush(stdout);}
	pfda_gen_e1_ab( zeta_ab, lambda, inv_sigma_eta, ka, kb, dl);
	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("zetaab:\n");printmat(zeta_ab,*ka,*kb);fflush(stdout);}
	
	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("compute zeta_bb:\n");fflush(stdout);}
	pfda_gen_e1_bb( zeta_bb, inv_sigma_eta, psi, xi, M, n, kb, dl);
	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("zeta_bb:\n");printmat(zeta_bb,*kb,*kb);fflush(stdout);}
	
	pfdaFree(inv_sigma_eta);
	if(checkdebug(dl,debugnum_dual_e1)){pfda_debug_msg("leaving pfdaDual_e1\n");fflush(stdout);}
}

void pfda_gen_e2_1(
	int    const * const outer,
	int    const * const a1,
	int    const * const c1,
	double const * const A,
	double const * const B,
	double       * const C,
	double       * const D,
	int const * const dl, double * dp,int * ip)
/* DESCRIPTION
INPUTS:
	a1
	c1
	A (a1 x a1)
	B (a1 x c1) or (c1 x a1)
	C (c1 x c1)  symetric		
OUTPUT:
	symetric matrices are returned in the upper portion
	D=(A-B C\inv C\trans)\inv      if outer 
	D=(A-B\trans C\inv \trans)    otherwise
	C=C\inv
MEMORY:
	dp	length = max( 2 a1^2, 2 c1^2)
	ip	length =  max ( a1, c1 )
dl: 281
*/
{//code
	if(checkdebug(dl,debugnum_pfda_gen_e2_1)){
		pfda_debug_msg("pfdaDual_e2_1 - \n");
		pfda_debug_msg("outer:\n%d=%s\n\n",*outer,(*outer)?("use outer"):("use inner"));
		pfda_debug_msg("a1:\n%d\n\n",*a1);
		pfda_debug_msg("c1:\n%d\n\n",*c1);
		pfda_debug_msg("A:\n"); printmat(A, *a1, *a1);
		pfda_debug_msg("B:\n"); if(outer) printmat(B, *a1, *c1); else printmat(B, *c1, *a1);
		pfda_debug_msg("C:\n"); printmat(C, *c1, *c1);
		fflush(stdout);
	}
	int sr;
	int a2=*a1**a1;

	//invert C
	double * C_inv = C;
	pfda_sym_inverse(C_inv, c1, &sr, dl, dp, ip);
	if(sr){pfda_error("PFDA ERR - pfdaDual_e2_1: Leading minor of order %i is not positive definite (C)",sr);}
	if(checkdebug(dl,debugnum_pfda_gen_e2_1)){pfda_debug_msg("C_inv:\n"); printmat(C_inv, *c1, *c1); fflush(stdout);}
		
	if(*outer){ pfda_matrix_outer_quadratic_form(D, B, a1, a1, C_inv, c1, dl, dp);
	} else {    pfda_matrix_inner_quadratic_form(D, B, a1, c1, C_inv, c1, dl, dp); }
	if(checkdebug(dl,debugnum_pfda_gen_e2_1)){
		if(*outer){pfda_debug_msg("B(C_inv)B^t:\n"); printmat(D, *a1, *a1); 
		}else{     pfda_debug_msg("B^t(C_inv)B:\n"); printmat(C_inv, *a1, *a1);}
		fflush(stdout);
	}
	dscal_(&a2, &mOne, D,&one); 
	daxpy_(&a2, &dOne, A, &one, D, &one);
	if(checkdebug(dl,debugnum_pfda_gen_e2_1)){pfda_debug_msg("A-i/o_q_form(B,C):\n"); printmat(D, *a1, *a1); fflush(stdout);}

	pfda_sym_inverse(D, a1, &sr, dl, dp, ip);
	if(sr){pfda_error("PFDA ERR - pfdaDual_e2_1: Leading minor of order %i is not positive definite (D)",sr);}
	if(checkdebug(dl,debugnum_pfda_gen_e2_1)){pfda_debug_msg("D=solve(A-i/o_q_form(B,C)):\n"); printmat(D, *a1, *a1); fflush(stdout);}
	if(checkdebug(dl,debugnum_pfda_gen_e2_1)){pfda_debug_msg("leaving pfda_gen_e2_1\n"); fflush(stdout);}
}

void pfda_gen_e2_aa(
	double       * const sigma_aa,
	double const * const zeta_aa,
	double const * const zeta_ab,
	double       * const zeta_bb,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,

	double * dp, int * ip)
/* DESCRIPTION
	(zeta_aa - zeta_ab (zeta_bb)\inv zeta_ab\trans)\inv
MEMORY
 	dp	length = max( 2 ka^2, 2 kb^2)
	ip	length =  max ( ka, kb )
*/
{
	pfda_gen_e2_1(&one,  ka, kb, zeta_aa, zeta_ab, zeta_bb, sigma_aa, dl, dp, ip);
}

void pfda_gen_e2_bb(
	double * const sigma_bb,
	double       * const zeta_aa,
	double const * const zeta_ab,
	double       * const zeta_bb,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,

	double * dp, int * ip)
/* DESCRIPTION
	(zeta_aa - zeta_ab (zeta_bb)\inv zeta_ab\trans)\inv
MEMORY
 	dp	length = max( 2 ka^2, 2 kb^2)
	ip	length =  max ( ka, kb )
*/
{
	pfda_gen_e2_1(&zero, kb, ka, zeta_bb, zeta_ab, zeta_aa, sigma_bb, dl, dp, ip);
}

void pfda_gen_e2_ab(
	double       * const sigma_ab,
	double const * const sigma_aa,
	double const * const zeta_ab,
	double const * const zeta_bb_inv,
	int    const * const ka,
	int    const * const kb,
	int const * const dl,
	double * dp)
/* DESCRIPTION
	sigma _ab = - sigma_aa zeta_ab (zeta_bb)inv
MEMORY
	dp length  = ka*kb
*/
{
	if(checkdebug(dl,debugnum_pfda_gen_e2_ab)){pfda_debug_msg("entering pfda_gen_e2_ab\n"); fflush(stdout);}
	double * tmp = pfdaAlloc_d(*ka**kb, &dp);
	dgemm_(&NoTrans, &NoTrans, ka,kb,ka, &dOne,  sigma_aa, ka, zeta_ab, ka, &dzero, tmp, ka);
	if(checkdebug(dl,debugnum_pfda_gen_e2_ab)){pfda_debug_msg("sigma_aa zeta_ab:\n"); printmat(tmp, *ka, *kb); fflush(stdout);}
	dsymm_(&Right, &Upper, ka, kb, &mOne, zeta_bb_inv, kb, tmp, ka, &dzero, sigma_ab, ka);	
	pfdaFree(tmp);
	if(checkdebug(dl,debugnum_pfda_gen_e2_ab)){pfda_debug_msg("leaving pfda_gen_e2_ab\n"); fflush(stdout);}
}

void pfdaDual_e2(
	int    const * const ka,
	int    const * const kb,
	double       * const zeta_aa,
	double const * const zeta_ab,
	double       * const zeta_bb,
	double       * const sigma_aa,
	double       * const sigma_ab,
	double       * const sigma_bb,
	int const * const dl, double * dp, int * ip	)
/* DESCRIPTION
INPUTS:
	ka
	kb
	zeta_aa (ka x ka)
	zeta_ab (ka x kb)
	zeta_bb (kb x kb)
OUTPUT:
	symetric matrices are returned in the upper portion
	\f[ %tex description
		\sigma_aa = (\zeta_aa-\zeta_ab \zeta_bb\inv \zeta_ab\trans)\inv
		\sigma_ab = -\sigma_aa\zeta_ab \zeta_bb\inv
		\sigma_bb = (\zeta_bb-\zeta_ab\trans \zeta_aa\inv \zeta_ab)\inv
	\f]
MEMORY:
	dp	length = kb^2 + max ( 2ka^2, 2*kb^2 )   use 3* max(ka,kb)^2
	ip	length = max ( ka ,kb)
*/
{//code
	if(checkdebug(dl,debugnum_dual_e2)){
		pfda_debug_msg("pfdaDual_e2 - \n");
		pfda_debug_msg("ka:\n%d\n\n",*ka);
		pfda_debug_msg("kb:\n%d\n\n",*kb);
		pfda_debug_msg("zeta_aa:\n"); printmat(zeta_aa, *ka, *ka);
		pfda_debug_msg("zeta_ab:\n"); printmat(zeta_ab, *ka, *kb);
		pfda_debug_msg("zeta_bb:\n"); printmat(zeta_bb, *kb, *kb);
		fflush(stdout);
	}
 	int kb2=*kb**kb;
	double * zeta_bb_inv = pfdaAlloc_d(kb2,&dp);
	
	dcopy_(&kb2,zeta_bb, &one, zeta_bb_inv, &one);
	if(checkdebug(dl,debugnum_dual_e2)){pfda_debug_msg("zeta_bb_inv:\n"); printmat(zeta_bb_inv, *kb, *kb); fflush(stdout);}
	pfda_gen_e2_aa( sigma_aa, zeta_aa, zeta_ab, zeta_bb_inv, ka, kb, dl, dp, ip);
	if(checkdebug(dl,debugnum_dual_e2)){pfda_debug_msg("zeta_bb_inv:\n"); printmat(zeta_bb_inv, *kb, *kb); fflush(stdout);}
	pfda_gen_e2_bb( sigma_bb, zeta_aa, zeta_ab, zeta_bb, ka, kb, dl, dp, ip);
	// pfdaDual_e2_1(&one,  ka, kb, zeta_aa, zeta_ab, zeta_bb_inv, sigma_aa, dl, dp, ip);
	// pfdaDual_e2_1(&zero, kb, ka, zeta_bb, zeta_ab, zeta_aa    , sigma_bb, dl, dp, ip);
	if(checkdebug(dl,debugnum_dual_e2)){pfda_debug_msg("sigma_aa:\n"); printmat(sigma_aa, *ka, *ka); fflush(stdout);}
	if(checkdebug(dl,debugnum_dual_e2)){pfda_debug_msg("sigma_bb:\n"); printmat(sigma_bb, *kb, *kb); fflush(stdout);}
	
	pfda_gen_e2_ab( sigma_ab, sigma_aa, zeta_ab, zeta_bb_inv, ka, kb, dl, dp);
	if(checkdebug(dl,debugnum_dual_e2)){
		pfda_debug_msg("sigma_ab:\n"); printmat(sigma_ab, *ka, *kb); 
		pfda_debug_msg("leaving pfdaDual_e2\n");fflush(stdout);}
}

void pfda_gen_e3_1(
	double       * const mu,
	int	   const * const ldmu,
	double const * const epsilon,
	double const * const AA,
	double const * const phi,
	int    const * const ldphi,
	double const * const u,
	int    const * const nu,
	int    const * const a,
	double const * const xi,
	double const * const AB,
	int    const * const transAB,
	double const * const psi,
	int    const * const ldpsi,
	double const * const v,
	int    const * const nv,
	int    const * const b,
	int const * const dl, double * dp)
	/* DESCRIPTION
	INPUTS:
		M- total number of observations
		epsilon  - singleton variance > 0
		AA (a x a) - symetric
		phi (nu x a) - with total number of rows N
		ldphi - the number of rows total for phi
		u,  - residuals
		nu - number of observations for u
		a,
		xi,
		AB (a x b) 
		transAB - boolean does AB have to be transposed?
		psi (nv x b), 
		ldpsi - the number of rows(total) for psi
		v,
		nv - number of observations for v
		b,
	OUTPUT:
		mu (a) = 1/epsilon * AA psi u + 1/xi * AB phi v
	MEMORY:
		dp	length = 2 a*n
	*/
{//code
	if(checkdebug(dl,debugnum_pfda_gen_e3_1)){
		pfda_debug_msg("entering pfda_gen_e3_1 - \n");
		pfda_debug_msg("epsilon:\n"); printmat(epsilon, one, one);
		pfda_debug_msg("AA:\n"); printmat(AA, *a, *a);
		pfda_debug_msg("phi:\n"); printmat2(*nu,*a, phi,ldphi);
		pfda_debug_msg("ldphi:\n%d\n\n",*ldphi);
		pfda_debug_msg("u:\n"); printmat(u, one, *nu);
		pfda_debug_msg("nu:\n%d\n\n",*nu);
		pfda_debug_msg("a:\n%d\n\n",*a);
		pfda_debug_msg("xi:\n"); printmat(xi, one, one);
		pfda_debug_msg("AB:\n"); printmat(AB, *a, *b);
		pfda_debug_msg("psi:\n"); printmat2(*nv,*b, psi,ldpsi);
		pfda_debug_msg("ldpsi:\n%d\n\n",*ldpsi);
		pfda_debug_msg("v:\n"); printmat(v, one, *nv);
		pfda_debug_msg("nv:\n%d\n\n",*nv);
		pfda_debug_msg("b:\n%d\n\n",*b);
		fflush(stdout);
	}
	double * Xa = pfdaAlloc_d(*a**nu,&dp);
	double * Xb = pfdaAlloc_d(*a**nv,&dp);
	double inv_epsilon = epsilon ? 1/(*epsilon) : 1;   if(checkdebug(dl,debugnum_pfda_gen_e3_1)){pfda_debug_msg("inv_epsilon:\n%9.8g\n\n",inv_epsilon);fflush(stdout);}
	double inv_xi      = xi      ? 1/(*xi)      : 1;             if(checkdebug(dl,debugnum_pfda_gen_e3_1)){pfda_debug_msg("inv_xi:\n%9.8g\n\n",inv_xi);fflush(stdout);}
	
	dgemm_(&NoTrans, &Trans, a, nu, a, &inv_epsilon, AA, a, phi, ldphi, &dzero, Xa, a);
	dgemm_((*transAB)?(&Trans):(&NoTrans), &Trans, a, nv, b, &inv_xi, AB, (*transAB)?b:a, psi, ldpsi, &dzero, Xb, a);
	if(checkdebug(dl,debugnum_pfda_gen_e3_1)){pfda_debug_msg("Xa:\n"); printmat(Xa, *a, *nu); fflush(stdout);}
	if(checkdebug(dl,debugnum_pfda_gen_e3_1)){pfda_debug_msg("Xb:\n"); printmat(Xb, *a, *nv); fflush(stdout);}
	
	dgemv_(&NoTrans, a, nu, &dOne, Xa, a, u, &one, &dzero, mu, ldmu); 
	if(checkdebug(dl,debugnum_pfda_gen_e3_1)){pfda_debug_msg("Xa u:\n"); printmat2(one, *a,mu,ldmu); fflush(stdout);}
	dgemv_(&NoTrans, a, nv, &dOne, Xb, a, v, &one, &dOne , mu, ldmu);
	pfdaFree(Xa);pfdaFree(Xb);
	if(checkdebug(dl,debugnum_pfda_gen_e3_1)){
		pfda_debug_msg("mu:\n"); printmat2(one, *a,mu,ldmu);
		pfda_debug_msg("leaving pfda_gen_e3_1\n");fflush(stdout);}
}

void pfda_dual_e3_1(
	double       * const mu,
	int	   const * const ldmu,
	int    const * const M,
	int    const * const n,
	double const * const epsilon,
	double const * const AA,
	double const * const phi,
	double const * const u,
	int    const * const a,
	double const * const xi,
	double const * const AB,
	int    const * const transAB,
	double const * const psi,
	double const * const v,
	int    const * const b,
	int const * const dl, double * dp)
	/* DESCRIPTION
	INPUTS:
		n - number of observations for subject
		M- total number of observations
		epsilon  - singleton variance > 0
		AA (a x a) - symetric
		phi (n x a) - with total number of rows N
		u, (a) - residuals
		a,
		xi,
		AB (a x b) 
		transAB - boolean does AB have to be transposed?
		psi (n x b), 
		v,
		b,
	OUTPUT:
		mu (a) = 1/epsilon * AA psi u + 1/xi * AB phi v
	MEMORY:
		dp	length = 2 a*n
	*/
{//code
	pfda_gen_e3_1( mu, ldmu, epsilon, AA, phi, M, u, n, a, xi, AB, transAB, psi, M, v,  n, b, dl, dp);
}

/*! Computes the variance matrices of alpha and beta given the data.

@MEMORY
	- dp = ka^2 + kb^2 + ka*kb + max(
		- kb2 + ka2 + max ( ka*kb  , 2*kb2)
		- kb^2 + max ( 2ka^2, 2*kb^2 )   use 3* max(ka,kb)^2
	7 * max(ka,kb)^2 is ample to cover it.
	- ip = max(ka, kb)
*/
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
	double*dp, int * ip)
{
	double * zeta_aa = pfdaAlloc_d(*ka**ka,&dp);
	double * zeta_ab = pfdaAlloc_d(*ka**kb,&dp);
	double * zeta_bb = pfdaAlloc_d(*kb**kb,&dp);
	pfda_dual_e1(M, ni, lambda, Da,  phi,   sep, ka, Db, psi , sxi, kb, zeta_aa, zeta_ab, zeta_bb, dl,dp, ip);
	pfdaDual_e2(ka,kb,zeta_aa,zeta_ab,zeta_bb,Saa,Sab,Sbb,dl,dp,ip);
}

void pfda_dual_e3(
	double       * const alpha,
	double       * const beta,
	double const * const sigma_aa,
	double const * const sigma_ab,
	double const * const sigma_bb,
	double const * const phi,
	double const * const psi,
	double const * const epsilon,
	double const * const xi,
	double const * const u,
	double const * const v,
	int    const * const ka,
	int    const * const kb,
	int    const * const n,
	int    const * const N,
	int    const * const M,
	int const * const dl, double * dp)
	/* DESCRIPTION
	INPUTS:
		sigma_aa,
		sigma_ab,
		sigma_bb,
		phi,
		psi,
		epsilon,
		xi,
		u,
		v,
		ka,
		kb,
		n - numer of obervations for this subject ( nrow or psi and phi)
		N - number of rows of matrix from which alpha  and beta are taken
		M - number of rows from which phi and psi are taken
	OUTPUT:
		mu (a) = 1/epsilon * sigma_aa phi u + 1/xi * sigma_ab phi v
		mu (a) = 1/epsilon * sigma_ab^t phi u + 1/xi * sigma_bb psi v
	MEMORY:
		dp	length= 2 *n * max(a,b)
	*/
{//code
	if(checkdebug(dl,debugnum_dual_e3)){pfda_debug_msg("Entering pfdaDual_e3\n");fflush(stdout);}
	pfda_dual_e3_1(alpha,N,M,n,epsilon,sigma_aa,phi,u,ka,xi     ,sigma_ab,&zero,psi,v,kb,dl,dp);
	pfda_dual_e3_1(beta ,N,M,n,xi     ,sigma_bb,psi,v,kb,epsilon,sigma_ab,&one ,phi,u,ka,dl,dp);
	if(checkdebug(dl,debugnum_dual_e3)){
		pfda_debug_msg("alpha:\n"); printmat2(one, *ka,alpha, N);
		pfda_debug_msg("beta: \n");  printmat2(one, *kb,beta, N);
		pfda_debug_msg("leaving pfdaDual_e3 \n");
		fflush(stdout);
	}
}

void pfdaDual_e(
	/* input Values */
	double const * const y, 
	double const * const z, 
	int    const * const nobs,
	int    const * const M, 
	int    const * const N, 
	int    const * const ka, 
	int    const * const kb, 
	const double * const B, 
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
	double       * const alpha, 
	double       * const Beta, 
	double       * const Sigma_aa, 
	double       * const Sigma_ab, 
	double       * const Sigma_bb,
	/* Control Values */
	int const * const dl,
	double * dp,int * ip)
/*DESCRIPTION:
	computes the Estep of the Paired fucntional data analysis problem to Zhou, et al.
MEMORY:
	dp	length = 2*M+M*ka+M*kb+ka^2+ka*kb+kb^2
		+PLUS+
		Max of:
			step 1.	kb^2 + ka^2 + max ( ka*kb  , 2*kb2)       < 4*max(ka,kb)^2
			step 2.	kb^2 + max ( 2ka^2, 2*kb^2 )                  < 3*max(ka,kb)^2
			step 3.	2 *max(ni) * max(a,b)
			= max( 2*max(ni)*max(a,b),4*max(ka,kb)^2)  might give a little extra, the but the extra should be really small provided ka and kb are small and/or close to each other.
	ip	length  = max ( ka ,kb)
*/
{ ///* E Step */
 	{ //debug
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("Entering pfdaDual_e\n");fflush(stdout);}
	if(checkdebug(dl,debugnum_dual_e)){
			pfda_debug_msg("INPUTS\n----------------------------------------\n");
			pfda_debug_msg("y:\n"); printmat(y, 1, *M);
			pfda_debug_msg("z:\n"); printmat(z, 1, *M);
			pfda_debug_msg("M:\n%d\n\n",*M);
			pfda_debug_msg("N:\n%d\n\n",*N );
			pfda_debug_msg("ka:\n%d\n\n",*ka);
			pfda_debug_msg("kb:\n%d\n\n",*kb);
			pfda_debug_msg("p:\n%d\n\n",*p);
			pfda_debug_msg("nobs:\n"); printmati(nobs, one, *N);
			pfda_debug_msg("tm:\n"); printmat(tm, 1, *p);
			pfda_debug_msg("tn:\n"); printmat(tn, 1, *p);
			pfda_debug_msg("tf:\n"); printmat(tf, *p, *ka);
			pfda_debug_msg("tg:\n"); printmat(tg, *p, *kb);
			pfda_debug_msg("Da:\n"); printmat(Da, 1, *ka);
			pfda_debug_msg("Db:\n"); printmat(Db, 1, *kb);
			pfda_debug_msg("Lambda:\n"); printmat(Lambda, *kb, *ka);
			pfda_debug_msg("seps:\n%g\n\n", *seps); 
			pfda_debug_msg("sxi:\n%g\n\n", *sxi); 
			pfda_debug_msg("alpha:\n"); printmat(alpha, *N, *ka);
			pfda_debug_msg("Beta:\n"); printmat(Beta, *N, *kb);
		}
	if(checkdebug(dl,debugnum_dual_e_inputs_large)){
			pfda_debug_msg("Sigma_aa:\n"); printmat(Sigma_aa, *N, *ka**ka);
			pfda_debug_msg("Sigma_aa:\n"); printmat(Sigma_aa, *N, *ka**ka);
			pfda_debug_msg("Sigma_ab:\n"); printmat(Sigma_ab, *N, *ka**kb);
			pfda_debug_msg("B:\n"); printmat(B, *M, *p);
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
	
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("Ry:\n"); printmat(Ry, one, *M); fflush(stdout);}
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("Rz:\n"); printmat(Rz, one, *M); fflush(stdout);}

	dgemm_(&NoTrans, &NoTrans, M, ka, p, &dOne, B, M, tf, p, &dzero, phi, M);//  phi = B*tf
	dgemm_(&NoTrans, &NoTrans, M, kb, p, &dOne, B, M, tg, p, &dzero, psi, M);//  psi = B*tg
	
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("phi:\n"); printmat(phi, *M, *ka); fflush(stdout);}
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("psi:\n"); printmat(psi, *M, *kb); fflush(stdout);}
		
	for(obsnum=0;obsnum<*N;obsnum++){
		if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("observation: %d ... ",obsnum+1);fflush(stdout);}
		if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("e1 ... "); fflush(stdout);}
		pfda_dual_e1(M,	nobs+obsnum, Lambda, 
				Da, phii, seps, ka,	
				Db, psii, sxi , kb, 
				zeta_aa, zeta_ab, zeta_bb,
				dl,dp,ip);
		if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("e2 ... "); fflush(stdout);}
		pfdaDual_e2(ka,kb,
				zeta_aa,zeta_ab,zeta_bb,
				Sigma_aa+obsnum**ka**ka,Sigma_ab+obsnum**ka**kb,Sigma_bb+obsnum**kb**kb,
				dl,dp,ip);
		if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("e3\n"); fflush(stdout);}
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
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("Freeing memory\n");fflush(stdout);}
	pfdaFree(Ry);
	pfdaFree(Rz);
	pfdaFree(phi);
	pfdaFree(psi);
	pfdaFree(zeta_aa);
	pfdaFree(zeta_ab);
	pfdaFree(zeta_bb);
	if(checkdebug(dl,debugnum_dual_e)){pfda_debug_msg("Leaving pfdaDual_e\n");fflush(stdout);}
	}
}


