unittestsdir = $(prefix)/unit-tests
unittests_SCRIPTS = run_unit_tests.sh

AM_FCFLAGS =                            \
        -I $(top_builddir)/src/         \
        -I $(top_builddir)/src/mpi	\
        -I $(top_builddir)/src/utils    \
        -I $(top_builddir)/src/fft      \
        -I $(top_builddir)/src/netcdf

AM_LDFLAGS =                                   	\
        $(top_builddir)/src/utils/*.o           \
        $(top_builddir)/src/fft/*.o          	\
        $(top_builddir)/src/fields/*.o          \
        $(top_builddir)/src/inversion/*.o

noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD =                            	\
	$(top_builddir)/src/mpi/libps_mpi.la       	\
        $(top_builddir)/src/utils/libps_utils.la      	\
        $(top_builddir)/src/fft/libps_mpi_fft.la       	\
        $(top_builddir)/src/netcdf/libps_netcdf.la


#         test_zeta
unittests_PROGRAMS =   		\
        test_vor2vel_1		\
        test_vor2vel_2		\
        test_vor2vel_3		\
        test_vor2vel_4		\
        test_vor2vel_5		\
        test_diffz_1		\
        test_diffz_2		\
        test_diffz_3		\
        test_diffz_4		\
        test_diffx		\
        test_diffy		\
        test_omp_collapse

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

test_diffz_1_SOURCES = unit_test.f90 test_diffz_1.f90
test_diffz_1_LDADD = libcombi.la

test_diffz_2_SOURCES = unit_test.f90 test_diffz_2.f90
test_diffz_2_LDADD = libcombi.la

test_diffz_3_SOURCES = unit_test.f90 test_diffz_3.f90
test_diffz_3_LDADD = libcombi.la

test_diffz_4_SOURCES = unit_test.f90 test_diffz_4.f90
test_diffz_4_LDADD = libcombi.la

test_diffx_SOURCES = unit_test.f90 test_diffx.f90
test_diffx_LDADD = libcombi.la

test_diffy_SOURCES = unit_test.f90 test_diffy.f90
test_diffy_LDADD = libcombi.la

# test_zeta_SOURCES = unit_test.f90 test_zeta.f90
# test_zeta_LDADD = libcombi.la

test_omp_collapse_SOURCES = unit_test.f90 test_omp_collapse.f90
test_omp_collapse_LDADD = libcombi.la
