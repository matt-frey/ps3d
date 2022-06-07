unittestsdir = $(prefix)/unit-tests
unittests_SCRIPTS = run_unit_tests.sh

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


unittests_PROGRAMS =   		\
        test_vor2vel_1		\
        test_vor2vel_2		\
        test_vor2vel_3		\
        test_vor2vel_4		\
        test_vor2vel_5

test_vor2vel_1_SOURCES = unit_test.f90 test_vor2vel_1.f90
test_vor2vel_1_LDADD = libcombi.la

test_vor2vel_2_SOURCES = unit_test.f90 test_vor2vel_2.f90
test_vor2vel_2_LDADD = libcombi.la

test_vor2vel_3_SOURCES = unit_test.f90 test_vor2vel_3.f90
test_vor2vel_3_LDADD = libcombi.la

test_vor2vel_4_SOURCES = unit_test.f90 test_vor2vel_4.f90
test_vor2vel_4_LDADD = libcombi.la

test_vor2vel_5_SOURCES = unit_test.f90 test_vor2vel_5.f90
test_vor2vel_5_LDADD = libcombi.la
