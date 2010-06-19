library(pfda)


#utilities
pfda:::UT_checkdebug()
pfda:::UT_pfdaAlloc()
pfda:::UT_pfda_eigens()
pfda:::UT_pfda_sym_inverse()
pfda:::UT_pfda_transpose()
pfda:::UT_pfda_uppersymadd()
pfda:::UT_pfda_computebtb()
pfda:::UT_pfda_matrix_inner_quadratic_form()
pfda:::UT_pfda_matrix_outer_quadratic_form()

#Random Generation
pfda:::UT_pfda_gen_truncnorm()

# binary utilities
pfda:::UT_pfda_bin_s_gen_w()
pfda:::UT_pfda_bin_single_approximate_moments_forobs()
pfda:::UT_pfda_bin_single_approximate_moments()
pfda:::UT_pfda_bin_cond_aa()
pfda:::UT_pfda_bin_cond_ab()
pfda:::UT_pfda_bin_cond_bb()

# single continuous steps
pfda:::UT_pfda_s_i()
pfda:::UT_pfda_m1()
pfda:::UT_pfda_m2()
pfda:::UT_pfda_m3()
pfda:::UT_pfda_m5_1()
pfda:::UT_pfdaSingle_e()

# single binary
pfda:::UT_pfda_bin_single_generate_w_parms1()
pfda:::UT_pfda_bin_single_e()

# dual functions
pfda:::UT_pfdaDual_e()
pfda:::UT_pfdaDual_e2()
pfda:::UT_pfdaDual_m4()
pfda:::UT_pfdaDual_m5_2()
pfda:::UT_pfda_cond_ab()
pfda:::UT_pfda_cond_dd()
pfda:::UT_pfda_dual_e1()
pfda:::UT_pfda_dual_e3_1()
pfda:::UT_pfda_dual_e3()
pfda:::UT_pfda_gen_e2_1()
pfda:::UT_pfda_sum_cond_dd()
pfda:::UT_dual_gen_sigmas()

# dual mixed bc
pfda:::UT_dual_bc_i()    # FAIL
pfda:::UT_dual_bc_1a()
pfda:::UT_dual_bc_1b()
pfda:::UT_dual_bc_1cde()
pfda:::UT_dual_bc_1()    # FAIL
pfda:::UT_dual_bc_2()    # FAIL
pfda:::UT_dual_bc_3()
pfda:::UT_dual_bc_4()
pfda:::UT_dual_bc_5()
pfda:::UT_dual_bc_6()
pfda:::UT_dual_bc_genw() 


eval(pfda:::X_dual_bc_genw)

