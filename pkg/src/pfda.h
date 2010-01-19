/*!  \file pfda.h main header for the pfda package compiled code.
 by Andrew Redd
//Part of paired functional data analysis through principle component analysis.
//This header contains definitions and constants used accross files.
*/
#ifndef PFDA_H
#define PFDA_H
#include <R.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <math.h>
#include "pfda_const.h"
#include "pfda_error.h"
#include "pfda_memory.h"
#include "pfda_utils.h"


#ifndef INLINE
# if __GNUC__
#  define INLINE __attribute__ ((gnu_inline))
# else
#  define INLINE inline
# endif
#endif

/// the only definitions allowed here are the \ref interfaces.

/*! \page interfaces Interface Functions
These functions are the interface functions.  The functions that R calls to enter the compiled portion of the code:
	- pfdaSingle()  the linear univariate
	- pfdaDual() the linear bivariate
	- pfda_bin_single() the binary univariate
*/	
#endif /*PFDA_H*/
