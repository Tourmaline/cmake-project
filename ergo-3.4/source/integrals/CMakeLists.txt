set(integrals_sources
basis_func_extent_1el.cc
basis_func_extent.cc
basis_func_pair_list_1el.cc
basis_func_pair_list.cc
box_system.cc
boysfunction.cc
densityfitting.cc
exponent_list.cc
#hermite_conversion.cc
hermite_conversion_explicit.cc
hermite_conversion_prep.cc
hermite_conversion_symb.cc
integral_info.cc
integrals_1el.cc
integrals_1el_kinetic.cc
integrals_1el_potential.cc
integrals_1el_single.cc
integrals_2el_boxed.cc
integrals_2el_coulomb.cc
integrals_2el_exchange.cc
integrals_2el_exchange_prep.cc
integrals_2el_explicit.cc
integrals_2el_layer.cc
integrals_2el_repeating.cc
integrals_2el_single.cc
integrals_2el_util_funcs.cc
integrals_2el_utils.cc
integrals_general.cc
integrals_hermite.cc
matrix_norm.cc
mm_limit_table.cc
monomial_info.cc
multipole.cc
operator_matrix.cc
organize_distrs.cc
simple_sparse_mat.cc
)

set(integrals_headers
basis_func_extent_1el.h
basis_func_extent.h
basis_func_pair_list_1el.h
basis_func_pair_list.h
box_system.h
boysfunction.h
densityfitting.h
exponent_list.h
hermite_conversion_explicit.h
#hermite_conversion.h
hermite_conversion_prep.h
hermite_conversion_symb.h
integral_info.h
integrals_1el.h
integrals_1el_kinetic.h
integrals_1el_potential.h
integrals_1el_single.h
integrals_2el_boxed.h
integrals_2el_coulomb.h
integrals_2el_exchange.h
integrals_2el_exchange_prep.h
integrals_2el_explicit.h
integrals_2el.h
integrals_2el_layer.h
integrals_2el_repeating.h
integrals_2el_single.h
integrals_2el_util_funcs.h
integrals_2el_utils.h
integrals_general.h
integrals_hermite.h
matrix_norm.h
mm_limit_table.h
monomial_info.h
multipole.h
operator_matrix.h
organize_distrs.h
simple_sparse_mat.h
)

add_library(integrals STATIC ${integrals_headers} ${integrals_sources})

# add subdirectories where to look for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/source)

set(library_subdirs
matrix
basisset
utilities
utilities_basic
)
 
foreach(dir ${library_subdirs})
	include_directories(${CMAKE_SOURCE_DIR}/source/${dir})
endforeach()

include_directories(${INTERNAL_BLAS_INCLUDE_DIRS})
include_directories(${INTERNAL_LAPACK_INCLUDE_DIRS})


install(TARGETS integrals
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${integrals_headers} DESTINATION include )



install(TARGETS integrals
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/static )
install( FILES ${integrals_headers} DESTINATION include )
