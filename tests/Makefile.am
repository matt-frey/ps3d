testsdir = $(prefix)/tests

AM_FCFLAGS =                            \
        -I $(top_builddir)/src/         \
        -I $(top_builddir)/src/mpi 	\
        -I $(top_builddir)/src/utils    \
        -I $(top_builddir)/src/fft      \
        -I $(top_builddir)/src/netcdf

AM_LDFLAGS =                                   	\
        $(top_builddir)/src/utils/*.o           \
        $(top_builddir)/src/fields/*.o          \
        $(top_builddir)/src/inversion/*.o

noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD =                            	\
	$(top_builddir)/src/mpi/libps_mpi.la            \
        $(top_builddir)/src/utils/libps_utils.la      	\
        $(top_builddir)/src/fft/libps_mpi_fft.la       	\
        $(top_builddir)/src/netcdf/libps_netcdf.la


tests_PROGRAMS =   		\
        test_vor2vel		\
        test_hyper
#         test_filter

test_vor2vel_SOURCES = test_vor2vel.f90
test_vor2vel_LDADD = libcombi.la

# test_filter_SOURCES = test_filter.f90
# test_filter_LDADD = libcombi.la

test_hyper_SOURCES = test_hyper.f90
test_hyper_LDADD = libcombi.la
