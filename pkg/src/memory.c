/*! \file memory.c

implimentation for the memory management functions
*/


#include "pfda.h"
///* New Memory allocation System */
int * pfdaAlloc_i(int size, int ** pool){
	int i;
	int * rtn = *pool;
	*pool+=size;
	for(i=0;i<size;i++)rtn[i]=0;
	return rtn;
}
double * pfdaAlloc_d(int size, double ** pool){
	int i;
	double * rtn = *pool;
	*pool+=size;
	for(i=0;i<size;i++)rtn[i]=0;
	return rtn;
}
void pfdaFree(void * p){}
void test_pfdaAlloc_int(int * PASS,int * dl){
	pfda_debug_cdl(debugnum_test_pfdaAlloc_int){pfda_debug_msg("test_pfdaAlloc_int\n");fflush(stdout);}
	int dp[100];
	int * dpp=dp;
	int * new;
	
	new = (int*)pfdaAlloc_i(10, &dpp);
	pfda_debug_cdl(debugnum_test_pfdaAlloc_int){
		pfda_debug_msg("dp:  %p\n", dp);
		pfda_debug_msg("dpp: %p\n", dpp);
		pfda_debug_msg("new: %p\n", new);
		fflush(stdout);
	}
	if(new==dp && dpp==dp+10) *PASS=1; else *PASS=0;	
}
void test_pfdaAlloc_double(int * PASS, int * dl){
	double dp[100];
	double * dpp=dp;
	double * new;
	
	new = (double*)pfdaAlloc_d(10, &dpp);
	pfda_debug_cdl(debugnum_test_pfdaAlloc_double){
		pfda_debug_msg("dp:  %p\n", dp);
		pfda_debug_msg("dpp: %p\n", dpp);
		pfda_debug_msg("new: %p\n", new);
		fflush(stdout);
	}
	if(new==dp && dpp==dp+10) *PASS=1; else *PASS=0;	
}


