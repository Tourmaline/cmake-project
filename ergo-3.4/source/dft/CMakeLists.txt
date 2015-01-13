set(dft_sources
aos.cc
cubature_rules.cc
dft_common.cc
fun-becke.c
fun-cam.c
functionals.c
fun-gga.c
fun-kt.c
fun-lb94.c
fun-lyp.c
fun-optx.c
fun-p86c.c
fun-pbec.c
fun-pbex.c
fun-pw86x.c
fun-pw91c.c
fun-pw91x.c
fun-pw92c.c
fun-pz81.c
fun-slater.c
fun-tester.cc
fun-vwn.c
grid_atomic.cc
grid_hicu.cc
grid_reader.cc
grid_stream.cc
integrator.cc
lebedev_laikov.cc
lin_trans.cc
rho-mat.cc
sparse_matrix.cc
sparse_pattern.cc
xc_matrix.cc
xc_matrix_sparse.cc
)

set(dft_headers
aos.h
cubature_rules.h
dft_common.h
dft.h
functionals.h
grid_atomic.h
grid_hicu.h
grid_interface.h
grid_matrix.h
grid_params.h
grid_reader.h
grid_stream.h
integrator.h
lebedev_laikov.h
lin_trans.h
rho-mat.h
sparse_matrix.h
sparse_pattern.h
xc_matrix.h
xc_matrix_sparse.h
)

add_library(dft STATIC ${dft_headers} ${dft_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
basisset
matrix
integrals
utilities
utilities_basic
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})



install(TARGETS dft
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${dft_headers} DESTINATION include )