# this line is needed such that config.h is found
AM_FCFLAGS = -I $(top_builddir)/src/mpi         \
	     -I $(top_builddir)/src/utils       \
             -I $(top_builddir)/src/netcdf

lib_LTLIBRARIES = libps_utils.la
libps_utils_la_SOURCES = 	\
	iomanip.f90		\
	constants.f90		\
	config.f90		\
	physics.f90		\
	jacobi.f90		\
	linsolve.f90		\
	merge_sort.f90		\
	dimensions.f90

if ENABLE_BUOYANCY
libps_utils_la_SOURCES += 	\
	ape_density.f90
endif
