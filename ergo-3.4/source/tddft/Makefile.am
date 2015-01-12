noinst_PROGRAMS = tddft_test

tddft_test_SOURCES = \
	tddft.cc	\
	tddft.h		\
	tddft_test.cc

tddft_test_LDADD = \
	$(top_builddir)/source/scf/libscf.a		\
	$(top_builddir)/source/dft/libdft.a 		\
	$(top_builddir)/source/densfromf/libdensfromf.a	\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
	$(top_builddir)/source/integrals/libintegrals.a			\
	$(top_builddir)/source/matrix/libmat.a		\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a           \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
	$(top_builddir)/source/basisset/libbasisset.a

# use matrix for gblas...
AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/dft		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas	\
        -I$(top_srcdir)/source/matrix/template_lapack/lapack	\
	-I$(top_srcdir)/source/scf		\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic

