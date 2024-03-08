SUBDIRS = mpi netcdf utils fft

AM_FCFLAGS =				\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/netcdf

bin_PROGRAMS = ps3d beltrami
#genspec genspec2d  production

ps3d_SOURCES = 					\
	utils/options.f90			\
	utils/parameters.f90			\
	utils/merge_sort.f90			\
	inversion/inversion_utils.f90		\
	fields/fields.f90			\
	fields/field_netcdf.f90			\
	fields/field_diagnostics.f90		\
	fields/field_diagnostics_netcdf.f90	\
	smagorinsky.f90				\
	inversion/inversion.f90			\
	utils/utils.f90				\
	advance.f90				\
	ps3d.f90

ps3d_LDADD = 						\
	$(top_builddir)/src/utils/libps_utils.la 	\
	$(top_builddir)/src/fft/libps_mpi_fft.la	\
	$(top_builddir)/src/mpi/libps_mpi.la		\
	$(top_builddir)/src/netcdf/libps_netcdf.la

#genspec_SOURCES = 				\
	utils/options.f90			\
	utils/parameters.f90			\
	inversion/inversion_utils.f90		\
	fields/fields.f90			\
	fields/field_netcdf.f90			\
	fields/field_diagnostics.f90		\
	fields/field_diagnostics_netcdf.f90	\
	inversion/inversion.f90			\
	utils/utils.f90				\
	genspec.f90

#genspec_LDADD = 					\
	$(top_builddir)/src/utils/libps_utils.la 	\
	$(top_builddir)/src/fft/libps_mpi_fft.la	\
	$(top_builddir)/src/mpi/libps_mpi.la		\
	$(top_builddir)/src/netcdf/libps_netcdf.la

#genspec2d_SOURCES = 				\
	utils/options.f90			\
	utils/parameters.f90			\
	inversion/inversion_utils.f90		\
	fields/fields.f90			\
	fields/field_netcdf.f90			\
	fields/field_diagnostics.f90		\
	fields/field_diagnostics_netcdf.f90	\
	inversion/inversion.f90			\
	utils/utils.f90				\
	genspec2d.f90

#genspec2d_LDADD = 					\
	$(top_builddir)/src/utils/libps_utils.la 	\
	$(top_builddir)/src/fft/libps_mpi_fft.la	\
	$(top_builddir)/src/mpi/libps_mpi.la		\
	$(top_builddir)/src/netcdf/libps_netcdf.la

beltrami_SOURCES = 				\
	utils/options.f90			\
	utils/parameters.f90			\
	inversion/inversion_utils.f90   	\
	fields/fields.f90			\
	fields/field_netcdf.f90			\
	fields/field_diagnostics.f90		\
	fields/field_diagnostics_netcdf.f90	\
        inversion/inversion.f90         	\
	utils/utils.f90				\
	beltrami.f90

beltrami_LDADD = 					\
	$(top_builddir)/src/utils/libps_utils.la 	\
	$(top_builddir)/src/fft/libps_mpi_fft.la        \
	$(top_builddir)/src/mpi/libps_mpi.la		\
	$(top_builddir)/src/netcdf/libps_netcdf.la

#production_SOURCES = 				\
	utils/options.f90			\
	utils/parameters.f90			\
	inversion/inversion_utils.f90   	\
	fields/fields.f90			\
	fields/field_netcdf.f90			\
	fields/field_diagnostics.f90		\
	fields/field_diagnostics_netcdf.f90	\
        inversion/inversion.f90         	\
	utils/utils.f90				\
	production.f90

#production_LDADD = 					\
	$(top_builddir)/src/utils/libps_utils.la 	\
	$(top_builddir)/src/fft/libps_mpi_fft.la	\
	$(top_builddir)/src/mpi/libps_mpi.la		\
	$(top_builddir)/src/netcdf/libps_netcdf.la

clean-local:
	rm -f *.mod
