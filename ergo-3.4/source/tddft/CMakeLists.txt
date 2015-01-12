set(tddft_sources
tddft.cc
tddft_test.cc
)

set(tddft_headers
tddft.h
)

add_executable(tddft_test ${tddft_headers} ${tddft_sources})
target_link_libraries(
    tddft_test
    scf
    densfromf
    dft
    utilities
    utilities_basic
    integrals
    mat
    basisset
    ${LAPACK_LIBRARIES}
    ${THREAD_LIBRARIES}
    )


# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
matrix
integrals
basisset
utilities
utilities_basic
scf
dft
densfromf
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})
