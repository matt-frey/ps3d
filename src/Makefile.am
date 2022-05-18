SUBDIRS = netcdf utils fft

AM_FCFLAGS =				\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/netcdf

bin_PROGRAMS = ps3d
ps3d_SOURCES = 				\
	utils/options.f90			\
	utils/parameters.f90			\
	fields/fields.f90			\
	fields/field_netcdf.f90			\
	inversion/inversion_utils.f90		\
	inversion/inversion.f90			\
	utils/utils.f90				\
	advance.f90				\
	ps3d.f90

ps3d_LDADD = 						\
	$(top_builddir)/src/utils/libps_utils.la 	\
	$(top_builddir)/src/fft/libps_fft.la		\
	$(top_builddir)/src/netcdf/libps_netcdf.la

clean-local:
	rm -f *.mod
