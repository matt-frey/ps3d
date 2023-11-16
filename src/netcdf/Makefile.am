AM_FCFLAGS = -I $(top_builddir)/src/ 		\
	     -I $(top_builddir)/src/utils 	\
	     -I $(top_builddir)/src/mpi

lib_LTLIBRARIES = libps_netcdf.la
libps_netcdf_la_SOURCES = 	\
	../utils/config.f90	\
	netcdf_utils.f90	\
	netcdf_writer.f90	\
	netcdf_reader.f90

libps_netcdf_la_LIBADD = $(top_builddir)/src/mpi/libps_mpi.la
