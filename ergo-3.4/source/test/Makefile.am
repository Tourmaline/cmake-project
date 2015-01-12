noinst_PROGRAMS = 		\
	basicmath_test		\
	lapack_test		\
	grid_test		\
	stream_test		\
	tmat_test		\
	vmat_test		\
	gradient_test		\
	mat_thread_test		\
	mat_performance_test	\
	jmat_test		\
	jmat_extrapolate_test	\
	kmat_test		\
	kmat_symm_test		\
	kmat_nosymm_test	\
	explicit_integral_test	\
	xcmat_test		\
	xcmat_sparse_test	\
	xcmat_nan_inf_test	\
	matmem_test		\
	mmul_simple_test	\
	mmul_simple_many_test	\
	mmul_rectangular_test

basicmath_test_SOURCES = basicmath_test.cc
lapack_test_SOURCES = lapack_test.cc
grid_test_SOURCES = grid_test.cc
stream_test_SOURCES = stream_test.cc
tmat_test_SOURCES = tmat_test.cc
vmat_test_SOURCES = vmat_test.cc
gradient_test_SOURCES = gradient_test.cc
mat_thread_test_SOURCES = mat_thread_test.cc
mat_performance_test_SOURCES = mat_performance_test.cc
jmat_test_SOURCES = jmat_test.cc
jmat_extrapolate_test_SOURCES = jmat_extrapolate_test.cc
kmat_test_SOURCES = kmat_test.cc
kmat_symm_test_SOURCES = kmat_symm_test.cc
kmat_nosymm_test_SOURCES = kmat_nosymm_test.cc
explicit_integral_test_SOURCES = explicit_integral_test.cc
xcmat_test_SOURCES = xcmat_test.cc
xcmat_sparse_test_SOURCES = xcmat_sparse_test.cc
xcmat_nan_inf_test_SOURCES = xcmat_nan_inf_test.cc
matmem_test_SOURCES = matmem_test.cc
mmul_simple_test_SOURCES = mmul_simple_test.cc
mmul_simple_many_test_SOURCES =  mmul_simple_many_test.cc
mmul_rectangular_test_SOURCES = mmul_rectangular_test.cc

basicmath_test_LDADD =	\
        $(top_builddir)/source/matrix/libmat.a          \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a

lapack_test_LDADD =	\
        $(top_builddir)/source/matrix/libmat.a          \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

grid_test_LDADD =	\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities/libutilities.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
        $(top_builddir)/source/matrix/libmat.a

stream_test_LDADD =	\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a \
        $(top_builddir)/source/matrix/libmat.a

tmat_test_LDADD =	\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a

vmat_test_LDADD =	\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a		\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

gradient_test_LDADD =	\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a		\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

mat_thread_test_LDADD =	\
	$(top_builddir)/source/utilities/libutilities.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

mat_performance_test_LDADD =	\
	$(top_builddir)/source/utilities/libutilities.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

jmat_test_LDADD =	\
	$(top_builddir)/source/scf/libscf.a			\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

jmat_extrapolate_test_LDADD =	\
	$(top_builddir)/source/scf/libscf.a			\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

kmat_test_LDADD =	\
	$(top_builddir)/source/scf/libscf.a			\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

kmat_symm_test_LDADD =	\
	$(top_builddir)/source/scf/libscf.a			\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

kmat_nosymm_test_LDADD =	\
	$(top_builddir)/source/scf/libscf.a			\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

explicit_integral_test_LDADD =	\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a

xcmat_test_LDADD =	\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
	$(top_builddir)/source/matrix/libmat.a

xcmat_sparse_test_LDADD =	\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
	$(top_builddir)/source/matrix/libmat.a

xcmat_nan_inf_test_LDADD =	\
	$(top_builddir)/source/dft/libdft.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
	$(top_builddir)/source/matrix/libmat.a

matmem_test_LDADD =	\
	$(top_builddir)/source/scf/libscf.a			\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

mmul_simple_test_LDADD =	\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

mmul_simple_many_test_LDADD =	\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

mmul_rectangular_test_LDADD =	\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a	\
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a	\
	$(top_builddir)/source/matrix/libmat.a

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset 	\
	-I$(top_srcdir)/source/dft		\
	-I$(top_srcdir)/source/scf		\
	-I$(top_srcdir)/source/densfromf	\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/matrix		\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack

TESTS = \
	basicmath_test	\
	lapack_test	\
	grid_test	\
	stream_test	\
	xcmat_test	\
	xcmat_sparse_test\
	xcmat_nan_inf_test	\
	tmat_test	\
	vmat_test	\
	gradient_test	\
	mat_thread_test	\
	mat_performance_test	\
	jmat_test	\
	jmat_extrapolate_test	\
	kmat_test		\
	kmat_symm_test		\
	kmat_nosymm_test	\
	explicit_integral_test	\
	matmem_test		\
	mmul_simple_test	\
	mmul_simple_many_test	\
	mmul_rectangular_test
