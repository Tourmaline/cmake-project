set(ci_sources
ci.cc	
)

set(ci_headers
ci.h
)

add_library(ci STATIC ${ci_headers} ${ci_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
basisset
densfromf
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

