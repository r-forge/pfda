/*! \file debug.c 
debugging functions 
*/
/*! \defgroup  debug Debugging
 *	Contains the functions for debugging.  
 *  The header file is pfda_error.h
 *  files should use the macro definitions which are turned on by defining the macrop variable DEBUG
 */
/*@{*/

#include "pfda.h"

// BEGIN FUNCTION DECLARATIONS //////////////////////////////////////////////////////////////
/*! check the debug level array dl agains the number level and returns a logical if debugging is active for that level.
 * the first element of dl should contain the length of the debug level array.
 * The array contains the length then indicators at each element.
 *  
 * the coller graph should be empty when debug is off since all functions should go through the macro pfda_debug_cdl
 *  @callergraph
 */
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
/*! Print a vector (double)
 *  Accessed through the macro pfda_debug_argvec
 *  @callergraph
 */
void printvec(double const * const y, int const * const n){
	for(int i=0; i<*n; i++)pfda_debug_msg("%# 9.8g\t", y[i]);
	pfda_debug_msg("\n");
}
/*! Print a vector (integer)
 *  Accesed through the macro pfda_debug_argveci
 * @callergraph
 */
void printveci(int const * const y, int const * const n){
	for(int i=0; i<*n; i++)pfda_debug_msg("%d\t", y[i]);
	pfda_debug_msg("\n");
}
/*! Print a vector of responses
 * each row represents a subject
  * accessed through macro pfda_debug_argyvec or through pfda_debug_yvec if nobs and N are the second and third arguments.
* @callergraph
 */
void printyvec(double const * const y, int const * const nobs, int const * const N){
	int i,j;
	double const * yi=y;
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
/*! print a vector of integer reponses
 * for use with binary model debugging
 * accessed through macro pfda_debug_argyveci or through pfda_debug_yveci if nobs and N are the second and third arguments.
 * @callergraph
 */
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
/*! Print a matrix 
 *  accessed through pfda_debug_argmat
 *	@callergraph
 */
void printmat(const double* M,const int nrow,const int ncol){//For Debug ONLY
int i,j;
for(i=0;i<nrow;i++){
	for(j=0;j<ncol;j++){
		pfda_debug_msg("%# 9.8g\t",M[j*nrow +i]);
	}
	pfda_debug_msg("\n");
}
// pfda_debug_msg("\n");
}
/*!	print a matrix of intergers
 *	accessed through pfda_debug_argmati
 *	@callergraph
 */
void printmati(const int const *M,const int nrow,const int ncol){//For Debug ONLY
int i,j;
for(i=0;i<nrow;i++){
	for(j=0;j<ncol;j++){
		pfda_debug_msg("%d\t",M[j*nrow +i]);
	}
	pfda_debug_msg("\n");
}
// pfda_debug_msg("\n");
}
/*! print a selection of a matrix
 * accessed through pfda_debug_argmat2
 * @callergraph
 */
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
/*! print a matrix storred in the packed mode.
 * @callergraph
 */
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
/*@}*/
