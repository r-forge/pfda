/*! \file pfdautils.c
Utility functions
*/
#include "pfda_utils.h"
#include "pfda_error.h"
#include "pfda_memory.h"
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifndef NULL
#define NULL 0
#endif

const char V='V';
const char Upper='U';
const char Lower='L';
const char NoTrans='N';
const char Trans = 'T';
const char NonUnit = 'N';
const char Unit = 'U';
const char Left = 'L';
const char Right = 'R';
const char NoVectors = 'N';
const char Vectors = 'V';
const char ValuesInInterval = 'V';
const char IthValues = 'I';
const char AllValues = 'A';
const int one=1;
const double dOne=  1.0;
const double mOne= -1.0;
const int zero=0;
const double dzero = 0.0;
const double MinVariance = 1E-16;

double* Idiag(int i){
	//pfda_debug_msg("Idiag is depreciated please use resetI with a previously allocated vector.\n");
	double* rtn;
	rtn = (double*)S_alloc(i*i,sizeof(double));
	int j;
	for(j=0;j<i;j++){
		rtn[j*i+j]=1.0;
	}
	return(rtn);
}
void resetI(double* A,const int n){
	int i,j;
	for(i=0;i<n;i++)for(j=0;j<n;j++){
		A[i*n+j]= (i==j) ? 1.0 : 0.0;
	}
	
}
int local_dspev(const int* n, double* A, double* W, double* Z, const int *ldz, double* work){
	//extern void dspev_(const char* jobz, const char* uplo, 
	//	const int* n, double * A, double* W, double* Z, const int* ldz, double* work, long* info);
	int info;
	dspev_(&V, &Upper, n, A, W, Z, ldz, work, &info);
	return info;
}
void positivefirstrow(double *A,const int nrow,const int ncol){
	int colnum,rownum;
	for(colnum=0;colnum<ncol;colnum++){
		if(A[colnum*nrow]<0)for(rownum=0;rownum<nrow;rownum++){
			A[nrow*colnum+rownum]*=-1.0;
		}
	}
}
void pfda_transpose(double * const A, const int nrow, const int ncol,
	int const * const dl, double * const dp)
{
	/*
	DESCRIPTION:
		Transposes matrix A.
	MEMORY:
		pd - nrow*ncol
	*/
	pfda_debug_cdl(debugnum_pfda_transpose){pfda_debug_msg("pfda_transpose - \n");fflush(stdout);}
	double * ldp = dp;
	if((ncol==1)||(nrow==1)) return;
	double * B=NULL; //temporary storage
	int i;
	int N=nrow*ncol;
	pfda_debug_cdl(debugnum_pfda_transpose){
		pfda_debug_msg("allocating B - \n");
		pfda_debug_arg(ldp);
		pfda_debug_arg(B);
		fflush(stdout);
	}
	B=pfdaAlloc_d(N,&ldp);
	pfda_debug_cdl(debugnum_pfda_transpose){
		pfda_debug_arg(ldp);
		pfda_debug_arg(B);
		pfda_debug_msg("B <- A\n");
		fflush(stdout);
	}
	dcopy_(&N,A, &one, B, &one);
	pfda_debug_cdl(debugnum_pfda_transpose){pfda_debug_msg("A <- t(B)\n");fflush(stdout);}
	for(i=0;i<nrow;i++){
		pfda_debug_cdl(debugnum_pfda_transpose){pfda_debug_msg("i: %d\n", i);fflush(stdout);}
		dcopy_(&ncol,B+i, &nrow, A+i*ncol, &one);
	}
	pfda_debug_cdl(debugnum_pfda_transpose){pfda_debug_msg("leaving pfda_transpose\n");fflush(stdout);}
}
void transpose(double *A, const int nrow, const int ncol)//depreciated
{
	/*
	DESCRIPTION:
		Transposes matrix A.
	MEMORY:
		pd - nrow*ncol
	*/
	//pfda_debug_msg("transpose is depreciated use pfda_transpose\n");
	if((ncol==1)||(nrow==1)) return;
	double * B; //temporary storage
	int i;
	int N=nrow*ncol;
	B=(double*)Calloc(N, double);
	dcopy_(&N,A, &one, B, &one);
	for(i=0;i<nrow;i++){
		dcopy_(&ncol,B+i, &nrow, A+i*ncol, &one);
	}
	Free(B);
}
void myinner(double* A,int* ncol,int* nrow,double* C,int* ldc){ //depreciated
	dsyrk_(&Upper, &Trans, ncol, nrow, &dOne, A, nrow, &dOne, C, ldc);
}
void _dqfmm(const char* UPLO, const int* m,const int* k,const double* alpha, const double* A, const int* lda, const double* S, const int* lds, const double* beta, double *Q, const int* ldq, const int debug){
	/*
	A is m by k
	S is k by k 
	q is m by m
	Q<-\alpha A S A' + \beta Q
	*/
	//pfda_debug_msg("_dqfmm is depreciated pleas use pfda_matrix_outer_quadratic_form.\n");
	if(debug){
		pfda_debug_msg("beta: %g\n", *beta);
		pfda_debug_argmat2(*m, *k, A, lda);
		pfda_debug_argmat2(*k, *k, S, lds);
		pfda_debug_argmat2(*m, *m, Q, ldq);
		fflush(stdout);
	}
	double* tmp=NULL;
	tmp = Calloc(*k**m, double);
	dsymm_(&Right, UPLO, m, k, &dOne, S, lds, A, lda, &dzero, tmp, m);//tmp=AS
	if(debug){
		pfda_debug_msg("tmp(A*S):\n");pfda_debug_argmat2(*m, *k, tmp, m);
	}
	dgemm_(&NoTrans, &Trans, m, m, k, alpha, tmp, m, A,lda, beta, Q, ldq);//Q=beta tmp A'=beta A S A'
	Free(tmp);
} 
void _dtqimm(const char* UPLO, const int* m,const int* k,const double* alpha, const double* A, const int* lda, double* S, const int* lds, const double* beta, double *Q, const int* ldq, int *info){
	/*
	A is m by k
	S is m by m
	Q is k by k
	Q<-\alpha A' S^-1 A + \beta Q
	On exit:
	Q contains the new quadratic form.
	S contains the cholesky decomposition. for S
	*/
	pfda_warning("_dqfmm is depreciated pleas use pfda_matrix_inner_quadratic_form.\n");
	int row, col;
	double* tmp=NULL;
	tmp = Calloc(*m**k, double); // m by k
	for(row=0;row<*k;row++)for(col=0;col<*m;col++)tmp[col+row**m]=A[col+row**lda];//copy A into tmp.
	dposv_( UPLO, m, k, S, lds, tmp, m, info ); //Tmp= S^-1 A
	if(!*info)dgemm_(&Trans, &NoTrans, k, k, m, alpha, A, lda, tmp, m, beta, Q, ldq);//Q=alpha tmp A'=alpha A S^-1 A'+beta Q
	Free(tmp);
} 
void _dqimm(const char* UPLO, const int* m,const int* k,const double* alpha, const double* A, const int* lda, double* S, const int* lds, const double* beta, double *Q, const int* ldq, int *info){
	/*
	A is m by k
	S is k by k
	Q is m by m
	Q<-\alpha A S^-1 A' + \beta Q
	On exit:
	Q contains the new quadratic form.
	S contains the cholesky decomposition. for S
	*/
	int i, j;
	double* tmp=NULL;
	tmp = Calloc(*m**k, double);// tmp is k by m
	for(i=0;i<*k;i++)for(j=0;j<*m;j++)tmp[j+i**m]=A[j+i**lda];//copy A into tmp.
	dposv_( UPLO, k, m, S, lds, tmp, k, info ); //Tmp= S^-1 A'
	if(!*info)dgemm_(&NoTrans, &NoTrans, m, m, k, alpha, A, lda, tmp, k, beta, Q, ldq);//Q=alpha A S^-1 A'+ beta Q
	Free(tmp);
} 

double tr(double * A,int n){
	double rtn = 0;
	int i;
	for(i=0;i<n;i++){
		rtn += A[i*(n)+i];
	}
	return rtn;
}
void setMinVar(double * const D, int const * const k, double const*const var){
	int i;
	for(i=0;i<*k;i++)D[i]=MMAX(D[i],*var);
}
void reverse(const int *N, double *V){
	int i;
	double tmp;
	for(i=0;i<*N/2;i++){
		tmp=V[i];
		V[i]=V[*N-i-1];
		V[*N-i-1]=tmp;
	}
}
void reverseColumns(const int *ncol, double *A, const int* lda){
	int i;
	for(i=0;i<*ncol/2;i++)dswap_(lda, A+i**lda, &one, A+*lda*(*ncol-i-1), &one);
	
}
void addDiag(int const * const N, double const * const D, double * const A, int const *const lda){
	/*
	D is an N vector describing a diagonal of a matrix.
	A is an lda by N matrix.
	the elements od D are added to the leading diagonal elements of A
	*/
	int element;
	for(element=0;element<*N;element++){
		A[element+element**lda]+=D[element];
	}
}
void addInvDiag(int const * const N, double const * const D, double * const A, int const *const lda){
	/*
	D is an N vector describing a diagonal of a matrix.
	A is an lda by N matrix.
	the elements od D are added to the leading diagonal elements of A
	*/
	int i;
	for(i=0;i<*N;i++){
		A[i+i**lda]+=1/D[i];
	}
}
void _dchsm(const char *side, const char *uplo, const int *M, const int *N,const double *alpha, const double *A, const int *lda, double *B, const int *ldb){
	/*
	solve from cholesky decomposition
	side - side to apply inverse
	uplo - upper or lower diagonal in A
	M - number of rows of B
	N - number of columns of B
	alpha - constant to multiply result by,
	A -  M by M matrix containing Cholesky decomposition in uplo triangle
	lda - first dimention of A
	B - matrix on side
	ldb - first dimention of B
	*/

	dtrsm_(side, uplo, ((*side==Left) ^ (*uplo==Upper))?&NoTrans:&Trans, &NonUnit, M, N, &dOne, A, lda, B, ldb);
	dtrsm_(side, uplo, ((*side==Left) ^ (*uplo==Upper))?&Trans:&NoTrans, &NonUnit, M, N, alpha, A, lda, B, ldb);
}
void uppersym2full(int const * const n, double const * const sym, double *const full, int const * const dl){
	pfda_debug_cdl(debugnum_uppersym2full){pfda_debug_msg("uppersym2full -");fflush(stdout);}
	int i,j;
	for(i=0;i<*n;i++)for(j=i;j<*n;j++){
		pfda_debug_cdl(debugnum_uppersym2full){pfda_debug_msg("sym[%d,%d]=%# 9g",i,j,sym[i]);fflush(stdout);}
		if(i==j) full[j**n+i]=sym[j**n+i]; else full[i**n+j]=full[j**n+i]=sym[j**n+i];
	}
	pfda_debug_cdl(debugnum_uppersym2full){pfda_debug_msg("leaving uppersym2full");fflush(stdout);}
}
void pfda_fillsym(double * const sym,int const * const n, int const * const dl){
	pfda_debug_cdl(debugmun_pfda_fillsym){pfda_debug_msg("uppersym2full -");fflush(stdout);}
	int i,j;
	for(i=0;i<*n;i++)for(j=i+1;j<*n;j++){
		pfda_debug_cdl(debugmun_pfda_fillsym){pfda_debug_msg("sym[%d,%d]=%# 9g",i,j,sym[i]);fflush(stdout);}
		sym[i**n+j]=sym[j**n+i];
	}
	pfda_debug_cdl(debugmun_pfda_fillsym){pfda_debug_msg("leaving uppersym2full");fflush(stdout);}
}
void test_LU_dchsm(const int *M, const int *N,const double *alpha, const double *A, const int *lda, double *B, const int *ldb){
	_dchsm(&Left, &Upper, M, N, alpha, A, lda, B, ldb);
}
void test_RU_dchsm(const int *M, const int *N,const double *alpha, const double *A, const int *lda, double *B, const int *ldb){
	_dchsm(&Right, &Upper, M, N, alpha, A, lda, B, ldb);
}
void test_LL_dchsm(const int *M, const int *N,const double *alpha, const double *A, const int *lda, double *B, const int *ldb){
	_dchsm(&Left, &Lower, M, N, alpha, A, lda, B, ldb);
}
void test_RL_dchsm(const int *M, const int *N,const double *alpha, const double *A, const int *lda, double *B, const int *ldb){
	_dchsm(&Right, &Lower, M, N, alpha, A, lda, B, ldb);
}

int pfda_arrMAXi(int * arr, int * length){
	int max=0;
	for(int i=0;i<*length;i++)
		max=MMAX(max,arr[i]);
	return max;
}


/*! \brief computes the outer quadratic form
 \f[
 Q = X S X^T
 \f]
	\param Q [output](nrowx x nrowx) symetric 
	\param X (nrowx x nrows)
	\param nrowx number of rows of X to use
	\param ldx the number of rows in the calling routine for X
	\param S (nrows by nrows) assumed symetric, only upper portion used.
	\param nrows number of rows of S
	\param dl
	\param dp see memory requirements
	dl: 306
	@MEMORY:
		dp - nrowx * nrows
*/
void pfda_matrix_outer_quadratic_form(
	double * const Q,
	double const * const X,
	int const * const nrowx,
	int const * const ldx,
	double const * const S,
	int const * const nrows,
	int const * const dl,
	double * dp )
{
	pfda_debug_cdl(debugnum_pfda_matrix_outer_quadratic_form){
		pfda_debug_msg("pfda_matrix_outer_quadratic_form - \n");
		pfda_debug_argmat(Q, *nrowx, *nrowx);
		pfda_debug_argmat(X, *nrowx, *nrows);
		pfda_debug_arg(*nrowx);
		pfda_debug_argmat(S, *nrows, *nrows);
		pfda_debug_arg(*nrows);
	}
	/* memory allocation */
	double* tmp=NULL;
	tmp = pfdaAlloc_d(*nrowx**nrows, &dp);
	
	/*!\f$tmp = XS\f$*/	dsymm_(&Right, &Upper, nrowx, nrows, &dOne, S, nrows, X, ldx, &dzero, tmp, nrowx);
	/*!\f$Q = tmp X^T = X S X^t\f$*/		dgemm_(&NoTrans, &Trans, nrowx, nrowx, nrows, &dOne, tmp, nrowx, X, ldx, &dzero, Q, nrowx);//Q=beta tmp A'=beta A S A'
	pfdaFree(tmp);
}

/*! computes the inner quadratic form
\f[	Q = X^T S X 	\f]

	\param	Q 	[output](ncolx x ncolx) symetric
	\param	X 	(nrows x ncols)
	\param	ncolx 	number of columns of X to use
	\param	ldx 	number of rows of X in the calling program to use
	\param	S 	(nrows by nrows) assumed symetric, only upper portion used.
	\param	nrows	the number of rows of S.
	\param	dl
	\param 	dp	 - ncolx*nrows
	dl: 307
	@MEMORY
	dp	 - ncolx*nrows
*/
void pfda_matrix_inner_quadratic_form(
		double       * const Q,
		double const * const X,
		int    const * const ncolx,
		int    const * const ldx,
		double const * const S,
		int    const * const nrows,
		int const * const dl,
		double * dp)
{
	pfda_debug_cdl(debugnum_pfda_matrix_inner_quadratic_form){
		pfda_debug_msg("pfda_matrix_inner_quadratic_form - \n");
		pfda_debug_argmat(Q, *ncolx, *ncolx);
		pfda_debug_argmat(X, *nrows, *ncolx);
		pfda_debug_arg(*ncolx);
		pfda_debug_argmat(S, *nrows, *nrows);
		pfda_debug_arg(*nrows);
	}
	///* memory allocation */
	double* tmp=NULL;
	tmp = pfdaAlloc_d(*ncolx**nrows, &dp);
	
	/*tmp = SX*/	dsymm_(&Left, &Upper, nrows, ncolx, &dOne, S, nrows, X, ldx, &dzero, tmp, nrows);
	pfda_debug_cdl(debugnum_pfda_matrix_inner_quadratic_form){pfda_debug_msg("tmp(S X):\n"); pfda_debug_argmat(Q, *nrows, *ncolx);}
	/*Q = X^T tmp = X^t S X*/		dgemm_(&Trans, &NoTrans, ncolx, ncolx, nrows, &dOne, X, ldx, tmp, nrows, &dzero, Q, ncolx	);//Q=beta tmp A'=beta A S A'
	pfdaFree(tmp);
	pfda_debug_cdl(debugnum_pfda_matrix_inner_quadratic_form){pfda_debug_argmat(Q, *ncolx, *ncolx);}
}


/*	returns the inverse of the symetric (upper) matrix A in A (a full matrix returned).
	 INPUTS:
		A ( n *n )
		n
	OUTPUTS:
		A = A^-1
		sr the result of the inversion.
	MEMORY:
		dp - 2 * n^2
		ip - n
	dl: 313
*/
void pfda_sym_inverse(double * const A,int const * const n, int * sr,
	int const * const dl, double * dp, int * ip)
{
	pfda_debug_cdl(debugnum_pfda_sym_inverse){
		pfda_debug_msg("pfda_sym_inverse -\n");
		pfda_debug_argmat(A,*n,*n);
		pfda_debug_msg("n: %d\n\n", *n);
		fflush(stdout); 
	}
	int n2=*n**n;
	double * inv = pfdaAlloc_d(n2,&dp);
	double * work = pfdaAlloc_d(n2,&dp);
	int * ipiv = pfdaAlloc_i(*n, &ip);
	resetI(inv,*n);
	dsysv_(&Upper, n,n, A, n, ipiv, inv, n, work, &n2, sr);
	// pfdaFree(ipiv);pfdaFree(work);
	dcopy_(&n2,inv,&one,A,&one);
	pfda_fillsym(A,n,dl);
	pfda_debug_cdl(debugnum_pfda_sym_inverse){
		pfda_debug_msg("A^-1:\n");pfda_debug_argmat(A,*n,*n);
		pfda_debug_msg("sr: %d\n\n", *sr);
		pfda_debug_msg("leaving pfda_sym_inverse\n");
		fflush(stdout); 
	}
}
void pfda_uppersymadd(
	int const * const n,
	double * const x, int const*const ldx,
	double const * const alpha,
	double const * const y, int const*const ldy,
	int const * const dl)
{
	int i,j;
	for(i=0;i<*n;i++)for(j=0;j<=i;j++)x[i**ldx+j]+=y[i**ldx+j]**alpha;
}


///* testing wrapper functions */
void test_transpose(double * A, int const* nrow, int const* ncol){transpose(A,*nrow,*ncol);}
void test_pfda_transpose(double * A, int const* nrow, int const* ncol, int const * const dl, double * dp){pfda_transpose(A,*nrow,*ncol,dl, dp);}
void test_positivefirstrow(double*A, int *nrow, int* ncol){positivefirstrow(A,*nrow,*ncol);}

