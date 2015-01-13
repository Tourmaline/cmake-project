add_subdirectory(purification)

set(densfromf_sources
densfromf_full.cc
densfromf_general.cc
densfromf_sparse.cc
densfromf_stochastic.cc
puri_mmul_tests.cc
simple_purification.cc
)

set(densfromf_headers
densfromf_full.h
densfromf_general.h
densfromf_sparse.h
densfromf_stochastic.h
puri_mmul_tests.h
simple_purification.h
)

add_library(densfromf STATIC ${densfromf_headers} ${densfromf_sources})
target_link_libraries(densfromf ${BLAS_LIBRARIES} ${LAPACK_LIBRARY})

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

include_directories(${CMAKE_SOURCE_DIR}/source/densfromf/purification)



install(TARGETS densfromf
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${densfromf_headers} DESTINATION include )
