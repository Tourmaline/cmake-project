set(basisset_sources
basisinfo.cc	
basisset.cc
)

set(basisset_headers
basisinfo.h	
basisset.h
)

add_library(basisset STATIC ${basisset_headers} ${basisset_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
matrix
integrals
utilities
utilities_basic
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()


install(TARGETS basisset
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${basisset_headers} DESTINATION include )
