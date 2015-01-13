# add all subdirectories
set(subdirs
matrix
integrals
basisset
electron_dynamics
utilities
utilities_basic
ci
scf
dft
densfromf
tddft
#test
)

foreach(dir ${subdirs})
	add_subdirectory(${dir})
endforeach()



# add subdirectories where to look for headers

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(library_subdirs
matrix
integrals
basisset
electron_dynamics 
utilities
utilities_basic
ci
scf
dft
densfromf
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})


# create ERGO executable

set(ergo_sources
ergo_scripted.cc
slr.cc	
)

set(ergo_headers
ergo_scripted.h
license.h
matrix_typedefs.h
slr.h
)


find_package(BISON)
find_package(FLEX)
BISON_TARGET(Parser ergo_input_processor.y ${CMAKE_CURRENT_BINARY_DIR}/ergo_input_processor.c)
FLEX_TARGET(Scanner ergo_input_parser.l ${CMAKE_CURRENT_BINARY_DIR}/ergo_input_parser.c)
ADD_FLEX_BISON_DEPENDENCY(Scanner Parser)

add_executable (ergo ${ergo_headers} ${ergo_sources} ${BISON_Parser_OUTPUTS} ${FLEX_Scanner_OUTPUTS})

target_link_libraries(
    ergo
    scf
    ci
    densfromf
    dft
    utilities
    utilities_basic
    integrals
    electron_dynamics
    mat
    basisset 
    ${THREAD_LIBRARIES}
    ${LAPACK_LIBRARIES}
    )


# cmake -DCMAKE_INSTALL_PREFIX=/path ..
# make install    (alternatively $ make -j<N> install)
install(TARGETS ergo
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${ergo_headers} DESTINATION include )



# additional executables

# fermi_contact
set(fermi_contact_sources
fermi_contact.cc	
)

set(fermi_contact_headers
)

add_executable (fermi_contact ${fermi_contact_headers} ${fermi_contact_sources})

target_link_libraries(
    fermi_contact
    basisset
    dft 
    utilities
    utilities_basic
    integrals
    mat
    ${LAPACK_LIBRARIES}
    )


# bin2m
set(bin2m_sources
bin2m.cc	
)

set(bin2m_headers
)

add_executable (bin2m ${bin2m_headers} ${bin2m_sources})

target_link_libraries(
    bin2m
    basisset 
    utilities
    utilities_basic
    integrals
    )


#inch_test
