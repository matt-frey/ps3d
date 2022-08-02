testsdir = $(prefix)/tests

AM_FCFLAGS =                            \
        -I $(top_builddir)/src/         \
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
        $(top_builddir)/src/utils/libps_utils.la      	\
        $(top_builddir)/src/fft/libps_fft.la          	\
        $(top_builddir)/src/netcdf/libps_netcdf.la


tests_PROGRAMS =   		\
        test_vor2vel

test_vor2vel_SOURCES = test_vor2vel.f90
test_vor2vel_LDADD = libcombi.la
