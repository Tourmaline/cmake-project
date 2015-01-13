set(utilities_sources
density_description_file.cc
densitymanager.cc
density_projection.cc
integral_matrix_wrappers.cc
matrix_utilities.cc
matrix_utilities_gao.cc
)

set(utilities_headers
density_description_file.h
densitymanager.h
density_projection.h
filematrix.h
integral_matrix_wrappers.h
mat_acc_extrapolate.h
matrix_utilities_gao.h
matrix_utilities.h
)

add_library(utilities STATIC ${utilities_headers} ${utilities_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
basisset
densfromf
matrix
integrals
utilities_basic
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})



install(TARGETS utilities
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${utilities_headers} DESTINATION include )

