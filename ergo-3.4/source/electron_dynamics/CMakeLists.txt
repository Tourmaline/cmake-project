set(electron_dynamics_sources
tdhf_dynamics.cc
)

set(electron_dynamics_headers
electron_dynamics.h	
tdhf_dynamics.h	
)

add_library(electron_dynamics STATIC ${electron_dynamics_headers} ${electron_dynamics_sources})

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
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${CMAKE_SOURCE_DIR}/source/matrix/template_lapack/blas)
include_directories(${CMAKE_SOURCE_DIR}/source/matrix/template_lapack/lapack)


install(TARGETS electron_dynamics
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${electron_dynamics_headers} DESTINATION include )
