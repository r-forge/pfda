/*! \file debug.c 
debugging functions 
*/

#include "pfda.h"

// BEGIN FUNCTION DECLARATIONS //////////////////////////////////////////////////////////////
int checkdebug(int const * const dl, const int level){
	/* Explanation of new debug system
		the new debug system opperates on a vector system. a boolean vector, with the first elelment the length of the vector. 
		The deeper the debugging the longer the arrary has to be, but since it is integer and debugging should be done relativly rarely. 
		This should not cause any memory problems.
		the existing system will only look at the first element and be ok.

	*/

	return (dl)&&(dl[0]>=level)&&(dl[level]);
}
void test_checkdebug(int const * const dl, int const * const level, int * rtn){*rtn = checkdebug(dl, *level);}
void printyvec(double const * const y, int const * const nobs, int const * const N){
	int i,j;
	double const * yi=y;
	pfda_debug_msg("obs  | y\n");
	for(i=0;i<*N;i++){
		pfda_debug_msg("%5d| ",i);
		for(j=0;j<nobs[i];j++){
			pfda_debug_msg("%# 9.8g\t", *yi);
			yi++;
		}
		pfda_debug_msg("\n");
	}
	pfda_debug_msg("\n");
}
void printyveci(int const * const y, int const * const nobs, int const * const N){
	int i,j;
	int const * yi=y;
	pfda_debug_msg("obs  | y\n");
	for(i=0;i<*N;i++){
		pfda_debug_msg("%5d| ",i);
		for(j=0;j<nobs[i];j++){
			pfda_debug_msg("%5d\t", *yi);
			yi++;
		}
		pfda_debug_msg("\n");
	}
	pfda_debug_msg("\n");
}
void printmat(const double* M,const int nrow,const int ncol){//For Debug ONLY
int i,j;
for(i=0;i<nrow;i++){
	for(j=0;j<ncol;j++){
		pfda_debug_msg("%# 9.8g\t",M[j*nrow +i]);
	}
	pfda_debug_msg("\n");
}
pfda_debug_msg("\n");
}
void printmati(const int const *M,const int nrow,const int ncol){//For Debug ONLY
int i,j;
for(i=0;i<nrow;i++){
	for(j=0;j<ncol;j++){
		pfda_debug_msg("%d\t",M[j*nrow +i]);
	}
	pfda_debug_msg("\n");
}
pfda_debug_msg("\n");
}
void printmat2(const int nrow,const int ncol, const double* A, const int* lda){//For Debug ONLY
int i,j;
for(i=0;i<nrow;i++){
	for(j=0;j<ncol;j++){
		pfda_debug_msg("%5.4f\t",A[j**lda +i]);
	}
	pfda_debug_msg("\n");
}
pfda_debug_msg("\n");
}
void printspmat(double* M, int p){
int i,j;
for(i=0;i<p;i++){
	for(j=0;j<i+1;j++){
		pfda_debug_msg("%5.4f\t",M[j*i+j]);
	}
	pfda_debug_msg("\n");
}
pfda_debug_msg("\n");
}
