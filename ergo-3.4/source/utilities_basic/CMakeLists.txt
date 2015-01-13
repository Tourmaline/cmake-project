set(utilities_basic_sources
atom_labels.c
barrier.c
csr_matrix.cc
ergo_math.cc
machine_epsilon.cc
matrix_algebra.cc
memorymanag.cc
molecule.cc
output.cc
solve_lin_eq_syst.cc
utilities.cc
xyz_file_parser.cc
)

set(utilities_basic_headers
atom_labels.h
barrier.h
csr_matrix.h
ergo_math.h
machine_epsilon.h
matrix_algebra.h
memorymanag.h
molecule.h
output.h
pi.h
polydegree.h
realtype.h
solve_lin_eq_syst.h
units.h
utilities.h
xyz_file_parser.h
)

add_library(utilities_basic STATIC ${utilities_basic_headers} ${utilities_basic_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
matrix
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})


install(TARGETS utilities_basic
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${utilities_basic_headers} DESTINATION include )
