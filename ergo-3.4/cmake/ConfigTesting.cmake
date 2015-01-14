
# activate ctest
include(CTest)
enable_testing()

# define a test
add_test(test_dft_hicu ${PROJECT_BINARY_DIR}/test_ctest/test_dft_hicu.sh)
add_test(test_6_d_funcs ${PROJECT_BINARY_DIR}/test_ctest/test_6_d_funcs.sh)
add_test(test_hf ${PROJECT_BINARY_DIR}/test_ctest/test_hf.sh)
add_test(test_ci ${PROJECT_BINARY_DIR}/test_ctest/test_ci.sh)
add_test(test_dft_hybrid ${PROJECT_BINARY_DIR}/test_ctest/test_dft_hybrid.sh)
add_test(test_dft_pure ${PROJECT_BINARY_DIR}/test_ctest/test_dft_pure.sh)
add_test(test_dft_sparse ${PROJECT_BINARY_DIR}/test_ctest/test_dft_sparse.sh)
add_test(test_dipole_moment ${PROJECT_BINARY_DIR}/test_ctest/test_dipole_moment.sh)
add_test(test_ext_elec_field ${PROJECT_BINARY_DIR}/test_ctest/test_ext_elec_field.sh)

