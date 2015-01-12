noinst_LIBRARIES = libutilities.a

libutilities_a_SOURCES = \
	integral_matrix_wrappers.cc	\
	integral_matrix_wrappers.h	\
	matrix_utilities.cc		\
	matrix_utilities.h		\
	matrix_utilities_gao.cc		\
	matrix_utilities_gao.h		\
	mat_acc_extrapolate.h		\
	density_description_file.cc	\
	density_description_file.h	\
	density_projection.cc		\
	density_projection.h		\
	densitymanager.cc		\
	densitymanager.h		\
	filematrix.h

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/densfromf	\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack
