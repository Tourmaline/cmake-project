SUBDIRS = \
	basisset 	\
	ci		\
	matrix 		\
	densfromf 	\
	dft 		\
	electron_dynamics	\
	integrals 	\
	scf 		\
	utilities	\
	utilities_basic	\
	.		\
	tddft		\
	test

bin_PROGRAMS = ergo
noinst_PROGRAMS = fermi_contact bin2m

AM_CPPFLAGS = \
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/dft 		\
	-I$(top_srcdir)/source/basisset 	\
	-I$(top_srcdir)/source/integrals 	\
	-I$(top_srcdir)/source/electron_dynamics 	\
	-I$(top_srcdir)/source/scf 		\
	-I$(top_srcdir)/source/ci 		\
	-I$(top_srcdir)/source/densfromf	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack

EXTRA_DIST = \
	$(top_srcdir)/source/standalone/standalone_density_file_conversion_program.c	\
	$(top_srcdir)/source/standalone/standalone_density_file_conversion_program_2.c

# ergo section
BUILT_SOURCES = ergo_input_processor.h
AM_YFLAGS = -d

ergo_SOURCES = \
	ergo_input_processor.y	\
	ergo_input_parser.l	\
	ergo_scripted.cc	\
	ergo_scripted.h		\
	license.h		\
	matrix_typedefs.h	\
	slr.cc			\
	slr.h

ergo_LDADD = \
	$(top_builddir)/source/scf/libscf.a		\
	$(top_builddir)/source/ci/libci.a		\
	$(top_builddir)/source/densfromf/libdensfromf.a	\
	$(top_builddir)/source/dft/libdft.a 		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a	\
	$(top_builddir)/source/integrals/libintegrals.a			\
	$(top_builddir)/source/electron_dynamics/libelectron_dynamics.a	\
	$(top_builddir)/source/matrix/libmat.a		\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a           \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
	$(top_builddir)/source/basisset/libbasisset.a

bin2m_SOURCES = bin2m.cc

bin2m_LDADD =	\
	$(top_builddir)/source/utilities/libutilities.a		\
	$(top_builddir)/source/basisset/libbasisset.a		\
	$(top_builddir)/source/integrals/libintegrals.a		\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a


fermi_contact_SOURCES = \
	fermi_contact.cc

fermi_contact_LDADD = \
	$(top_builddir)/source/dft/libdft.a 		\
	$(top_builddir)/source/utilities/libutilities.a	\
	$(top_builddir)/source/matrix/libmat.a		\
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a           \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a       \
	$(top_builddir)/source/basisset/libbasisset.a	\
	$(top_builddir)/source/integrals/libintegrals.a	\
	$(top_builddir)/source/utilities_basic/libutilities_basic.a


#ergo_input_parser.o : ergo_input_parser.c ergo_input_processor.h
