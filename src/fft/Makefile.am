AM_FCFLAGS = -I $(top_builddir)/src/ 		\
	     -I $(top_builddir)/src/utils 	\
	     -I $(top_builddir)/src/mpi

lib_LTLIBRARIES = libps_fft.la
libps_fft_la_SOURCES = 		\
	../utils/constants.f90	\
	deriv1d.f90		\
	stafft.f90		\
	sta2dfft.f90

lib_LTLIBRARIES += libps_mpi_fft.la

libps_mpi_fft_la_SOURCES = 		\
	../utils/constants.f90		\
	stafft.f90			\
	sta2dfft.f90			\
	fft_pencil.f90			\
	deriv1d.f90			\
	sta3dfft.f90

libps_mpi_fft_la_LIBADD = $(top_builddir)/src/mpi/libps_mpi.la 		\
			  $(top_builddir)/src/utils/libps_utils.la
