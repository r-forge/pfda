/*! 	\file eigens.c
	\brief Eigen Decomposition Functions
	
	Depends on LAPACK 
*/
#include "pfda_error.h"
#include "pfda_memory.h"
#include "pfda_utils.h"
#include <R_ext/Lapack.h>
#include "pfda_const.h"

/*!	Solves the eigen value problem:
	\ingroup utils
		\f[
			S => Q S Q^T
		\f]
		@param S (ns x ns) Square Symetrix matrix metrix, upper portion used. Destroyed on exit
		@param ns - number of rows and columns of S
		@param vectors [output] (ns x k) 
		@param values [output] (k)
		@param k  - number of eigenvalues to obtain
		@param dl - integer array for printing debugging information.
		@param dp see memory
		@param ip see memory
	@MEMORY:
		- dp - 8*ns
		- ip - 6*ns
*/
void pfda_eigens(
	double       * const S,  //destroyed on run
	int    const * const ns,
	double       * const vectors,	double       * const values,
	int    const * const k,
	int    const * const dl, double * dp, int * ip)
{

	pfda_debug_cdl(debugnum_pfda_eigens){
		pfda_debug_msg("pfda_eigens - \n");
		pfda_debug_msg("S:\n"); printmat(S, *ns, *ns);
		pfda_debug_msg("ns:\n%d\n\n",*ns );
		pfda_debug_msg("k:\n%d\n\n",*k );
		fflush(stdout);
	}
	int e_result=0;
	int lower_val_num=*ns-*k+1;
	int N_val_found=0;
	int lwork=8**ns;
	double * work=pfdaAlloc_d(lwork, &dp);	
	int * iwork=pfdaAlloc_i(5**ns,&ip);	
	int * ifail=pfdaAlloc_i(*ns, &ip);
	pfda_debug_cdl(debugnum_pfda_eigens){		
		pfda_debug_msg("lower_val_num:\n%d\n\n", lower_val_num );
		pfda_debug_msg("entering dsyevx:\n");
		fflush(stdout);
	}
	dsyevx_(&Vectors, &IthValues, &Upper, 
		ns, S, ns, NULL, NULL, &lower_val_num, ns, &dzero, &N_val_found, values, vectors, ns,
		work, &lwork, iwork, ifail, &e_result);
	pfda_debug_cdl(debugnum_pfda_eigens){		pfda_debug_msg("returned from dsyevx\n");fflush(stdout);}
	if(e_result)pfda_error("PFDA ERR: Could not find eigen values associated with D_alpha, info:%d",e_result);
	reverse(k, values);
	reverseColumns(k, vectors, ns);
	positivefirstrow(vectors, *ns, *k);	
	pfdaFree(work);
	pfdaFree(iwork);
	pfdaFree(ifail);
	pfda_debug_cdl(debugnum_pfda_eigens){
		pfda_debug_msg("vectors:\n"); printmat(vectors, *ns, *k);
		pfda_debug_msg("values:\n"); printmat(values, one, *k);
		pfda_debug_msg("Leaving pfda_eigens\n");
		fflush(stdout);
	}
}

/*! Computes the orthogonalization to find Da and tranforms tf and alpha

@MEMORY
	- dp = 9*p + 3*p^2 + N*ka
	- ip = 6p
*/
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
	double * dp, int * ip)
{
	double * S = pfdaAlloc_d(*p**p,&dp);
	pfda_matrix_outer_quadratic_form(S,tf,p,p,avgaa,ka,dl,dp );
	
	double * vectors = pfdaAlloc_d(*p**p,&dp);
	double * values  = pfdaAlloc_d(*p   ,&dp);
	pfda_eigens( S,p,vectors, values, ka, dl, dp, ip);
	positivefirstrow(vectors, *p, *ka);
	
	dcopy_(ka,values,&one,Da,&one);
	
	dgemm_(&Trans,&NoTrans, ka, ka, p, &dOne, vectors, p, tf, p, &dzero, trans, ka);
	
	int nka = *N**ka;
	double * new_alpha =pfdaAlloc_d(nka,&dp);
	dgemm_(&NoTrans,&Trans, N, ka, ka, &dOne, alpha  , N, trans, ka, &dzero, new_alpha, N);
	dcopy_(&nka, new_alpha, &one, alpha, &one);	
	
	int pka = *p**ka;
	dcopy_(&pka, vectors, &one, tf, &one);
}




