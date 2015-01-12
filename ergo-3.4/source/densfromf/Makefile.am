SUBDIRS = \
        purification

noinst_LIBRARIES = libdensfromf.a

libdensfromf_a_SOURCES = \
	densfromf_full.cc	\
	densfromf_full.h	\
	densfromf_general.cc	\
	densfromf_general.h	\
	densfromf_sparse.cc	\
	densfromf_sparse.h	\
	densfromf_stochastic.cc	\
	densfromf_stochastic.h	\
	puri_mmul_tests.h	\
	puri_mmul_tests.cc	\
	simple_purification.h	\
	simple_purification.cc

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas           \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack \
	-I$(top_srcdir)/source/densfromf/purification
