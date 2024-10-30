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
        test_omp_collapse	\
        test_zops		\
        test_clenshaw		\
        test_clenshaw_1		\
        test_clenshaw_2		\
        test_clenshaw_3		\
        test_integrate		\
        test_zfilter		\
        test_cheb_poly		\
        test_crank_nicolson	\
        test_crank_nicolson_1	\
        test_crank_nicolson_2	\
        test_crank_nicolson_3	\
	test_implicit_rk

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

test_zops_SOURCES = unit_test.f90 test_zops.f90
test_zops_LDADD = libcombi.la

test_clenshaw_SOURCES = unit_test.f90 test_clenshaw.f90
test_clenshaw_LDADD = libcombi.la

test_clenshaw_1_SOURCES = unit_test.f90 test_clenshaw_1.f90
test_clenshaw_1_LDADD = libcombi.la

test_clenshaw_2_SOURCES = unit_test.f90 test_clenshaw_2.f90
test_clenshaw_2_LDADD = libcombi.la

test_clenshaw_3_SOURCES = unit_test.f90 test_clenshaw_3.f90
test_clenshaw_3_LDADD = libcombi.la

test_integrate_SOURCES = unit_test.f90 test_integrate.f90
test_integrate_LDADD = libcombi.la

test_zfilter_SOURCES = unit_test.f90 test_zfilter.f90
test_zfilter_LDADD = libcombi.la

test_cheb_poly_SOURCES = unit_test.f90 test_cheb_poly.f90
test_cheb_poly_LDADD = libcombi.la

test_crank_nicolson_SOURCES = unit_test.f90 test_crank_nicolson.f90
test_crank_nicolson_LDADD = libcombi.la

test_crank_nicolson_1_SOURCES = unit_test.f90 test_crank_nicolson_1.f90
test_crank_nicolson_1_LDADD = libcombi.la

test_crank_nicolson_2_SOURCES = unit_test.f90 test_crank_nicolson_2.f90
test_crank_nicolson_2_LDADD = libcombi.la

test_crank_nicolson_3_SOURCES = unit_test.f90 test_crank_nicolson_3.f90
test_crank_nicolson_3_LDADD = libcombi.la

test_implicit_rk_SOURCES = unit_test.f90 test_implicit_rk.f90
test_implicit_rk_LDADD = libcombi.la
