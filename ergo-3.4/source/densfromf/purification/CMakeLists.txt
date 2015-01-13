set(test_sources
test.cc
)

set(test_headers
Purification_scaled.h	
Step.h		
TestMatrix.h	
)

add_executable(test ${test_headers} ${test_sources})


target_link_libraries(
    test
    mat
    templateblas
    templatelapack
    )

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
matrix
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${CMAKE_SOURCE_DIR}/source/matrix/template_lapack/blas)
include_directories(${CMAKE_SOURCE_DIR}/source/matrix/template_lapack/lapack)


