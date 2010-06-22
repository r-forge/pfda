/*! \file truncnorm.c

part of pfda package by Andrew Redd
Portions of this code are taken from 
Heike Trautmann, Detlef Steuer and Olaf Mersmann (2008). truncnorm:
  Trunctated Normal Distribution. R package version 1.0.0.
  http://r-forge.r-project.org/projects/desire
*/

#include "pfda_utils.h"
#include "pfda_error.h"
#include "pfda_memory.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <math.h>
#include <R_ext/Random.h>

void pfda_rtruncnorm(
	double       * const x,
	int    const * const n,
	double const * const lower, int    const * const nl,
	double const * const upper,	int    const * const nu,
	double const * const mean ,	int    const * const nm,
	double const * const sd   ,	int    const * const ns,
	int const * const dl)
{
	pfda_debug_cdl(debugnum_pfda_rtruncnorm){
		pfda_debug_msg("entering pfda_rtruncnorm\n");
		pfda_debug_msg("nl:%p\n",nl);
		pfda_debug_msg("nu:%p\n",nu);
		fflush(stdout);
		if(lower)pfda_debug_msg("lower: %g\n", *lower);
		if(upper)pfda_debug_msg("upper: %g\n", *upper);
		if(mean)pfda_debug_msg("mean: %g\n", *mean);
		if(sd)pfda_debug_msg("sd: %g\n", *sd);
		fflush(stdout);
	}
	if(nl && nu && *upper - *lower < dzero )pfda_error("error in pfda_rtruncnrom. Bad endpoints specified.");
	GetRNGstate();
	for (int i = 0; i < *n; ++i) {
		pfda_debug_cdl(debugnum_pfda_rtruncnorm){pfda_debug_msg("attempting to generate %d\n", i);fflush(stdout);}
		const double ca = nl?lower[i % *nl]:0.0 ;
		const double cb = nu?upper[i % *nu]:0.0 ;
		const double cmean = mean[i % *nm];
		const double csd = sd[i % *ns];
		double tmp;
		/* Rejection sampling until we generate an observation that is not
		 * truncated. This can be slow for narrow and extreme intervals
		 * [a,b].
		 */
		const int maxattempts = 100;
		int nattempts = maxattempts;
		while (--nattempts) { 
			tmp = rnorm(cmean, csd);
			pfda_debug_cdl(debugnum_pfda_rtruncnorm_attempts){
				pfda_debug_msg("candidate[%d]: %g\n",maxattempts-nattempts ,tmp);fflush(stdout);}
			if ((!nl || ca >= tmp) && (!nu || tmp <= cb))break;
		}
		if(!nattempts)pfda_error("Error:pfda_rtruncnorm: Could not find a suitable value in %d attempts.",maxattempts);
		x[i] = tmp;
	}
	PutRNGstate();
	pfda_debug_cdl(debugnum_pfda_rtruncnorm){pfda_debug_msg("leaving pfda_rtruncnorm\n");fflush(stdout);}
}
void pfda_rtruncnorm_from_below(
	double       * const x,
	int    const * const n,
	double const * const lower, int    const * const nl,
	double const * const mean ,	int    const * const nm,
	double const * const sd   ,	int    const * const ns,
	int const * const dl){ pfda_rtruncnorm(x, n, lower, nl, NULL, NULL, mean , nm, sd, ns, dl);}
void pfda_rtruncnorm_from_above(
	double       * const x,
	int    const * const n,
	double const * const upper,	int    const * const nu,
	double const * const mean ,	int    const * const nm,
	double const * const sd   ,	int    const * const ns,
	int const * const dl){ pfda_rtruncnorm(x, n, NULL, NULL, upper,nu, mean, nm, sd, ns, dl);}
	
/*! Approximates the truncated from below normal  distribution

@MEMORY
	* dp 	6 * m
	* ip	m
*/
void pfda_robert_truncnorm(
	double       * const x,
	int    const * const m,
	double const * const mean,
	double const * const sd,
	double const * const lower,
	int    const * const niterations,
	int const * const dl,
	double * dp,
	int * ip
)
{
	pfda_debug_cdl(debugnum_pfda_robert_truncnorm){pfda_debug_msg("niterations: %d\n",*niterations);fflush(stdout);}
	double * ct = pfdaAlloc_d(*m,&dp);
	double * R = pfdaAlloc_d(*m,&dp);
	double * U = pfdaAlloc_d(*m,&dp);
	double * Z = pfdaAlloc_d(*m,&dp);
	double * T1 = pfdaAlloc_d(*m,&dp);
	double * T2 = pfdaAlloc_d(*m,&dp);
	int    * s = pfdaAlloc_i(*m,&ip);
	int    * M = pfdaAlloc_i(*m,&ip);
	
	GetRNGstate();
	int i=0;
	double const 		c = (*lower-*mean)/(*sd);
	for(i=0;i<*m;i++){	ct[i]=x[i]=c;}
	double         		a = (c+sqrt(c*c+4));
	for(i=0;i<*m;i++){	s[i]  = ct>0;
						x[i] *= (double)s[i];
						R[i] = rnorm(0.0,1.0);
	}
	for(int j=0;j<*niterations;j++){
		pfda_debug_cdl(debugnum_pfda_robert_truncnorm){pfda_debug_msg("iteration: %d\n",j);fflush(stdout);}
		for(i=0;i<*m;i++)	U[i]= runif(0.0,1.0);
		for(i=0;i<*m;i++)	Z[i]= ct[i]-log(U[i])/a;
		for(i=0;i<*m;i++)	U[i]= runif(0.0,1.0);
		for(i=0;i<*m;i++)	M[i]= U[i] < exp(-(Z[i]-a)*(Z[i]-a)/2.0);
		for(i=0;i<*m;i++)	T1[i]= Z[i]*(double)M[i]+x[i]*(!(double)M[i]);
		for(i=0;i<*m;i++)	T2[i]= rnorm(0.0,1.0);
		for(i=0;i<*m;i++)	R[i] = R[i]*(T2[i]<=ct[i])+T2[i]*(T2[i]>ct[i]);
		for(i=0;i<*m;i++)	x[i] = R[i]*(ct[i])+T1[i]*(ct[i]>0);
	}
	PutRNGstate();
}

void pfda_test_random(
	double * const x,
	int const * const n)
{
	GetRNGstate();
	for(int i=0;i<*n;i++)x[i]=rnorm(0.0,1.0);
	PutRNGstate();
}

void pfda_roberts_1(
	double * const x,
	double const * const c,
	int const * const maxiterations,
	int const * const dl)
{
	//debugnum_pfda_roberts_1 = 431
	double const zero=0.0;
	double const one=1.0; 
	double a = (*c+sqrt(*c**c+4))/2;
	
	double p,v;
	int I = *maxiterations;
	GetRNGstate();
	do {
		double u = unif_rand();//runif(zero,one);
		*x = *c-log(u)/a;
		p = exp(-(*x-a)*(*x-a)/2.0) ;
		v = runif(zero,one);
		pfda_debug_cdl(debugnum_pfda_rtruncnorm_attempts){pfda_debug_msg("Attempts remaining = %d, x=%g, p=%g, v=%g\n",I,*x,p,v);}
		if(!(--I))pfda_error("pfda_robets_1: Excceeded maximum number of attempts to find a suitible random number.");
	} while(v>p);
	PutRNGstate();
}
void pfda_rejection_1(
	double * const x,
	double const * const c,
	int const * const max,
	int const * const dl
)
{
	int I=*max;
	GetRNGstate();
	do{
		*x = norm_rand();
		if(!(--I))pfda_error("Failed to find acceptable value with alloted number of attempts\n");
	} while(*x<*c);
	PutRNGstate();
}

/*! generates the tuncated normal using either the roberts or the rejection sampler depending on the tuncations value c. 
*/
void pfda_gen_truncnorm(
	double * x,            ///< [out]  the vector of generated values
	int const * const n,   ///< lenth of x to be generated
	double const * const c,///< the truncation value.  Truncation is from the left
	int const * const dl)  ///< the debug level
{
	pfda_debug_cdl(debugnum_pfda_gen_truncnorm){
		pfda_debug_msg("entering pfda_gen_truncnorm\n");
		pfda_debug_msg("n:%d\n",*n);
		pfda_debug_msg("c:%g\n",*c);
	}
	const int maxI=50;///the maximum number of itterations is build in and specified at 50.
	if(*c<0){
		for(int i=0;i<*n;i++)pfda_rejection_1(x+i, c, &maxI,dl);
	} else {
		for(int i=0;i<*n;i++)pfda_roberts_1(x+i, c, &maxI,dl);	
	}
}

void test_pfda_gen_truncnorm(
	double * x,
	int const * const n,
	double const * const c,
	int const * const dl)
{
GetRNGstate();
	pfda_gen_truncnorm(x,n,c,dl);
PutRNGstate();
}


