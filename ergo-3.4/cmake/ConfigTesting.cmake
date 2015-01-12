
# activate ctest
include(CTest)
enable_testing()

# define a test
add_test(test_dft_hicu ${PROJECT_BINARY_DIR}/test/test_dft_hicu.sh)
add_test(test_hf ${PROJECT_BINARY_DIR}/test/test_hf.sh)
