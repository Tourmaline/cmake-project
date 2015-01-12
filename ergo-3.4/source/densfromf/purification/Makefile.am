noinst_PROGRAMS = test

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack

test_SOURCES = \
	Purification_scaled.h	\
	Step.h		\
	TestMatrix.h	\
	test.cc	

test_LDADD   =      \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a \
        $(top_builddir)/source/matrix/libmat.a

