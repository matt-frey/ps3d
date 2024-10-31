# this line is needed such that config.h is found
AM_FCFLAGS = -I $(top_builddir)/src/mpi         \
	     -I $(top_builddir)/src/utils       \
             -I $(top_builddir)/src/netcdf	\
             -I $(top_builddir)/src/fft

lib_LTLIBRARIES = libps_utils.la
libps_utils_la_SOURCES = 	\
	iomanip.f90		\
	constants.f90		\
	parameters.f90		\
	config.f90		\
	physics.f90		\
	jacobi.f90		\
	merge_sort.f90		\
	dimensions.f90		\
	rolling_mean.f90	\
	cheby.f90

if ENABLE_BUOYANCY
libps_utils_la_SOURCES += 	\
	ape_density.f90
endif
