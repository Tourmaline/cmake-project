set(scf_sources
diis_general.cc
diis_restricted.cc
diis_unrestricted.cc
scf.cc
SCF_general.cc
SCF_restricted.cc
SCF_statistics.cc
SCF_unrestricted.cc
scf_utils.cc
)

set(scf_headers
diis_general.h
diis_restricted.h
diis_unrestricted.h
SCF_general.h
scf.h
SCF_restricted.h
SCF_statistics.h
SCF_unrestricted.h
scf_utils.h
)

add_library(scf STATIC ${scf_headers} ${scf_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
basisset
dft
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



install(TARGETS scf
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${scf_headers} DESTINATION include )