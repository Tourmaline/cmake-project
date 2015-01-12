set(templateblas_source
template_blas_basicmath.cc
template_blas_common.cc
)

set(templateblas_headers
template_blas_asum.h
template_blas_axpy.h
template_blas_basicmath.h
template_blas_common.h
template_blas_copy.h
template_blas_dot.h
template_blas_gemm.h
template_blas_gemv.h
template_blas_ger.h
template_blas_idamax.h
template_blas_nrm2.h
template_blas_rot.h
template_blas_scal.h
template_blas_spmv.h
template_blas_spr2.h
template_blas_spr.h
template_blas_swap.h
template_blas_symm.h
template_blas_symv.h
template_blas_syr2.h
template_blas_syr2k.h
template_blas_syrk.h
template_blas_tpmv.h
template_blas_tpsv.h
template_blas_trmm.h
template_blas_trmv.h
template_blas_trsm.h
template_blas_trsv.h
)

add_library(templateblas STATIC ${templateblas_headers} ${templateblas_source})


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source/)



add_executable(templateblastest template_blas_test.cc)
target_link_libraries(
templateblastest
templateblas
)


add_executable(templateblastest_threaded template_blas_test_threaded.cc)
target_link_libraries(
templateblastest_threaded
templateblas
${THREAD_LIBRARIES}
)
