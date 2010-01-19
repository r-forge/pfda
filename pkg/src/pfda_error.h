/*! 	@file pfda_error.h
	@brief debugging definitions, macros and debugnum enumeration
*/
#ifndef PFDA_ERROR_H_
#define PFDA_ERROR_H_
#include <R.h>

/*! \defgroup error Error and Debugging functions*/
/*\@{*/


#define pfda_error Rf_error /// For error messages.  Raises an R generic error.  Works like printf
#define pfda_warning Rf_warning /// For warning messages.  Raises an R generic Warning.  works like printf.
#define pfda_debug_msg Rprintf /// for debugging messages.  works like cat.  Requires an fflush(stdout); to print imediatly.


/*!  	@enum pfda_debug_num
	@brief the enumerated debugging numbers
	
	The debug numbers provide a constant to supply to the checkdebug function that checks against the debugging level supplied by R. 
	
	This is part of the pfda debugging system.
*/
enum pfda_debug_num{
	debugnum_s_i_break_alpha = 97,
	debugnum_s_1_break_tm = 98,
	debugnum_s_i = 99,
	// pfdaSingle (100-199)
	// ---------------------------
	debugnum_singe_steps				= 100,
	debugnum_single_inputs 				= 101,
	debugnum_single_inputs_large 			= 102,
	debugnum_single_estimates			= 103,
	debugnum_single_convergence			= 108,
	debugnum_single_convergence_components 	= 109,
	debugnum_m1 						= 110,
	debugnum_m1_inputs 					= 111,
	debugnum_m1_inloop	 				= 112,
	debugnum_m2 						= 120,
	debugnum_m2_inputs_large 			= 121,
	debugnum_m3 						= 134,
	debugnum_m3_level1					= 135,
	debugnum_m3_level2					= 136,
	debugnum_m3_level3					= 137,
	debugnum_m5 						= 150,
	debugnum_m5_0 						= 151,
	debugnum_m5_1 						= 152,
	debugnum_single_e 					= 160,
	debugnum_single_e_inloop 		 	= 161,
	debugnum_single_finish 				= 190,
	// pfdaDual (200-299)
	// --------------------------
	debugnum_dual_steps 				= 200,
	debugnum_dual_inputs 				= 201,
	debugnum_dual_inputs_large 			= 202,
	debugnum_dual 						= 203,
	debugnum_dual_estimates 			     = 204,
	debugnum_dual_criteria 				= 208,
	debugnum_dual_criteria_components 	     = 209,
	debugnum_dual_m4 					= 240,
	debugnum_dual_m5 					= 250,
	debugnum_dual_m5_1 					= 251,
	debugnum_dual_m5_2 					= 252,
	debugnum_dual_e 					= 260,
	debugnum_dual_e_inputs_large  		= 261,
	debugnum_dual_e1 					= 270,
	debugnum_dual_e2 					= 280,
	debugnum_dual_e3 					= 290,
	debugnum_dual_e3_1 					= 291,
	// others --------
	debugnum_cond_dd 					= 301,
	debugnum_cond_ab 					= 302,
	debugnum_sum_cond_dd				= 303,
	debugnum_sum_cond_ab 				= 304,
	debugnum_pfda_eigens 					= 305,
	debugnum_pfda_matrix_outer_quadratic_form 	= 306,
	debugnum_pfda_matrix_inner_quadratic_form 	= 307,
	debugnum_uppersym2full 				= 308 ,
	debugmun_pfda_fillsym                   = 309,
	debugnum_pfda_transpose 				= 310,
	debugnum_test_pfdaAlloc_int 			= 311,
	debugnum_test_pfdaAlloc_double 		= 312,
	debugnum_pfda_sym_inverse 			= 313,
	debugnum_pfda_computebtb 			= 314,
	debugnum_pfda_uppersymadd 			= 315,
	debugnum_pfda_rtruncnorm                = 316,
	debugnum_pfda_rtruncnorm_attempts       = 317,
	debugnum_pfda_matrix_drop1              = 318,
	// binary conditional,
	debugnum_pfda_bin_cond_aa               = 355,
	debugnum_test_pdfa_bin_cond_aa 		= 356,
	debugnum_pfda_bin_cond_bb 			= 357,
	debugnum_test_pfda_bin_cond_bb 		= 358,
	debugnum_pfda_bin_cond_ab 			= 359,
	// Bin and General
	debugnum_memory 					= 400,
	debugnum_pfda_gen_e_eta = 401,
	debugnum_pfda_gen_e1_aa = 402,
	debugnum_pfda_gen_e1_ab = 403,
	debugnum_pfda_gen_e1_bb = 404,
	debugnum_pfda_gen_e2_1  = 405,
	debugnum_pfda_gen_e2_aa = 406,
	debugnum_pfda_gen_e2_ab = 407,
	debugnum_pfda_gen_e2_bb = 408,
	debugnum_pfda_gen_e3_1  = 409,
	debugnum_pfda_bin_single_generate_w_parms1          = 411,
	debugnum_pfda_bin_single_generate_w_parms2          = 411,
	debugnum_pfda_bin_single_approximate_moments        = 412,
	debugnum_pfda_bin_single_approximate_moments_while  = 413,
	debugnum_pfda_bin_single_approximate_moments_forobs = 414,
	debugnum_bin_single_steps				     = 415,
	debugnum_bin_single_inputs 				= 416,
	debugnum_bin_single_inputs_large 			= 417,
	debugnum_bin_single_estimates			     = 418,
	debugnum_bin_single_convergence			= 419,
	debugnum_bin_single_convergence_components 	= 420,
	debugnum_pfda_bin_single_generate_w_1 	     = 421,
	debugnum_pfda_bin_m5 = 423,
	debugnum_pfda_bin_single_e = 424,
	debugnum_pfda_bin_single_break_before_e      = 425,
	debugnum_pfda_bin_single_break_before_2      = 426,
	debugnum_pfda_bin_single_break_before_3      = 427,
	debugnum_pfda_bin_single_break_before_4      = 428,
	debugnum_pfda_bin_single_break_before_ew     = 429,
	debugnum_pfda_bin_single_break_after_ew      = 430,
	debugnum_pfda_robert_truncnorm = 431,
	debugnum_pfda_roberts_1 = 432,
	debugnum_pfda_bin_s_gen_w = 433,
	// Calcium
	debugnum_dual_ca_inputs = 450
};
extern int  checkdebug(int const * const dl, const int level);
extern void printyvec(double const * const y, int const * const nobs, int const * const N);
extern void printyveci(int const * const y, int const * const nobs, int const * const N);
extern void printmat(const double* M,const int nrow, const int ncol);
extern void printmati(const int const *M,const int nrow, const int ncol);
extern void printmat2(const int nrow,const int ncol, const double* A, const int* lda);
extern void printspmat(double* M, int p);
/*\@}*/
#endif /* PFDA_ERROR_H_ */

