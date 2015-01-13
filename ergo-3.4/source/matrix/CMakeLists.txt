add_subdirectory(template_lapack)

set(mat_sources
allocate.cc
AllocatorManager.cc
FileWritable.cc
mat_gblas.cc
matInclude.cc
Memory_buffer_thread.cc
SizesAndBlocks.cc
#TCP_bench.cc
)

set(mat_headers
allocate.h
Allocator.h
AllocatorManager.h
bisection.h
DebugPolicies.h
Failure.h
FileWritable.h
general.h
Interval.h
Lanczos.h
LanczosLargestMagnitudeEig.h
mat_gblas.h
matInclude.h
MatrixBase.h
MatrixGeneral.h
Matrix.h
MatrixHierarchicBase.h
matrix_proxy.h
MatrixSymmetric.h
MatrixTriangular.h
MatrixTridiagSymmetric.h
mat_utils.h
Memory_buffer_thread.h
Perturbation.h
Purification.h
purification_old.h
PuriInfo.h
PuriStepInfoDebug.h
PuriStepInfo.h
SizesAndBlocks.h
sort.h
TC2.h
truncation.h
ValidPtr.h
VectorGeneral.h
Vector.h
VectorHierarchicBase.h
gemm_sse/common.h       
gemm_sse/g_intrin.h     
gemm_sse/mm_kernel_inner_sse2_A.h 
gemm_sse/mm_kernel_outer_A.h 
gemm_sse/gemm_sse.h 
gemm_sse/vector_intrin.h
)

add_library(mat STATIC ${mat_headers} ${mat_sources})
target_link_libraries(mat ${LAPACK_LIBRARY})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/gemm_sse)
include_directories(${CMAKE_SOURCE_DIR}/source/)

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})


# executables
add_executable (mattest API_test.cc)
target_link_libraries(
    mattest
    mat
    ${LAPACK_LIBRARIES}
    )

add_executable (matbench bench.cc)
target_link_libraries(
    matbench
    mat
    ${LAPACK_LIBRARIES}
    )

add_executable (matbenchgemm bench_gemm_only.cc)
target_link_libraries(
    matbenchgemm
    mat
    ${LAPACK_LIBRARIES}
    )

add_executable (omptest testOmpFail.cc)


add_executable (blastime blastime.cc)
target_link_libraries(
    blastime
    mat
    ${LAPACK_LIBRARIES}
    )

add_executable (perturbtest Perturb_Test.cc)
target_link_libraries(
    perturbtest
    mat
    ${LAPACK_LIBRARIES}
    )
	


install(TARGETS mat
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${mat_headers} DESTINATION include )




