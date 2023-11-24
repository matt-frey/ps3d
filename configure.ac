AC_INIT([ps3d], [0.1.1], [mf248@st-andrews.ac.uk], [], [https://github.com/matt-frey/ps3d])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_FC([gfortran])
AC_LANG(Fortran)
AM_PROG_AR
LT_INIT

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)

FCFLAGS="-std=f2018 -fdefault-real-8 -fdefault-double-8 -cpp -mcmodel=large"

AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([src/utils/config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    src/mpi/Makefile
    src/utils/Makefile
    src/netcdf/Makefile
    src/fft/Makefile
    unit-tests/Makefile
    tests/Makefile
])

#######################################################################################
##
## "--with" flags
##


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NetCDF C

## 11 March 2021
## https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
AC_ARG_WITH([netcdf-c],
            [AS_HELP_STRING([--with-netcdf-c=<path>], [path to NetCDF C library @<:@default=$NETCDF_C_DIR@:>@])],
            [NETCDF_C_DIR=$withval],
            [])

if test -n "${NETCDF_C_DIR}"; then
    CFLAGS="$CFLAGS -I$NETCDF_C_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_C_DIR/lib"
else
    AC_MSG_ERROR([No NETCDF_C_DIR environment variable.])
fi

AC_LANG_PUSH([C])
AC_SEARCH_LIBS([nc_open], [netcdf], [], [], [])
AC_MSG_CHECKING([whether we can compile a NetCDF C program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([], [
    #include <netcdf.h>
    int main() {
        int ncid;
        nc_create("dummy.nc", NC_CLOBBER, &ncid);
    }])],
[netcdf_c_found=yes],
[netcdf_c_found=no])
AC_MSG_RESULT([$netcdf_c_found])
AC_LANG_POP([C])

if test "x$netcdf_c_found" = "xno"; then
    AC_MSG_ERROR([Cannot compile a NetCDF C program])
fi


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NetCDF Fortran

AC_ARG_WITH([netcdf-fortran],
            [AS_HELP_STRING([--with-netcdf-fortran=<path>],
            [path to NetCDF Fortran library @<:@default=$NETCDF_FORTRAN_DIR@:>@])],
            [NETCDF_FORTRAN_DIR=$withval],
            [])

if test -n "${NETCDF_FORTRAN_DIR}"; then
    FCFLAGS="$FCFLAGS -I$NETCDF_FORTRAN_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_FORTRAN_DIR/lib"
else
    AC_MSG_ERROR([No NETCDF_FORTRAN_DIR environment variable.])
fi


# 15 April 2023
#http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_prog_fc_mpi.m4
# We cannot use AC_SEARCH_LIBS
for lib in NONE netcdff; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([for function nf90_open])
    else
    AC_MSG_CHECKING([for function nf90_open in -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [
           use netcdf
            implicit none
            integer :: ncerr, ncid
            ncerr = nf90_open("dummy.nc", NF90_NOWRITE, ncid)])],
        [netcdf_found=yes],
        [netcdf_found=no])
    AC_MSG_RESULT($netcdf_found)
    if test "x$netcdf_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done

AC_MSG_CHECKING([whether we can compile a NetCDF Fortran program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([], [
    use netcdf
    implicit none])],
[netcdf_fortran_found=yes],
[netcdf_fortran_found=no])
AC_MSG_RESULT([$netcdf_fortran_found])

if test "x$netcdf_fortran_found" = "xno"; then
    AC_MSG_ERROR([Cannot compile a NetCDF Fortran program])
fi


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MPI

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi=<path>], [path to mpi library @<:@default=$MPI_DIR@:>@])],
            [MPI_DIR=$withval],
            [])

if test -n "${MPI_DIR}"; then
    # mod file may also be in $MPI_DIR/lib
    FCFLAGS="$FCFLAGS -I$MPI_DIR/include -I$MPI_DIR/lib"
    LDFLAGS="$LDFLAGS -L$MPI_DIR/lib"
else
    AC_MSG_ERROR([No MPI_DIR environment variable.])
fi

AC_SEARCH_LIBS([MPI_Init],
               [mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90],
               [],
               [])

AC_SEARCH_LIBS([MPI_Win_lock],
               [mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90],
               [],
               [])

for lib in NONE mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([for module mpi_f08])
    else
    AC_MSG_CHECKING([for module mpi_f08 in -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
        [mpi_found=yes],
        [mpi_found=no])
    AC_MSG_RESULT([$mpi_found])

    if test "x$mpi_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done


AC_MSG_CHECKING([for MPI Fortran library])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
    [mpi_found=yes],
    [mpi_found=no])
AC_MSG_RESULT([$mpi_found])

if test "x$mpi_found" = "xno"; then
    AC_MSG_ERROR([Cannot find MPI Fortran library])
fi

#######################################################################################
##
## "--enable" flags
##

ENABLE_VERBOSE='no'
AC_ARG_ENABLE([verbose],
              [AS_HELP_STRING([--enable-verbose], [enable verbose mode (default=no)])],
              [ENABLE_VERBOSE=$enableval])

AM_CONDITIONAL([ENABLE_VERBOSE], [test "$ENABLE_VERBOSE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in verbose mode])
if test "x$ENABLE_VERBOSE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_VERBOSE"
else
    AC_MSG_RESULT([no])
fi

ENABLE_BUOYANCY='no'
AC_ARG_ENABLE([buoyancy],
              [AS_HELP_STRING([--enable-buoyancy], [enable buoyancy mode (default=no)])],
              [ENABLE_BUOYANCY=$enableval])

AM_CONDITIONAL([ENABLE_BUOYANCY], [test "$ENABLE_BUOYANCY" = "yes"])

AC_MSG_CHECKING([whether we are enabling buoyancy])
if test "x$ENABLE_BUOYANCY" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_BUOYANCY"
else
    AC_MSG_RESULT([no])
fi

ENABLE_PERTURBATION_MODE='no'
AC_ARG_ENABLE([perturbation-mode],
	      [AS_HELP_STRING([--perturbation-mode], [enable buoyancy perturbation mode (default=no)])],
              [ENABLE_PERTURBATION_MODE=$enableval])
AM_CONDITIONAL([ENABLE_PERTURBATION_MODE], [test "$ENABLE_PERTURBATION_MODE" = "yes"])
AC_MSG_CHECKING([whether we are tracking buoyancy perturbation])

if test "x$ENABLE_PERTURBATION_MODE" = "xyes"; then
   AC_MSG_RESULT([yes])
   FCFLAGS="$FCFLAGS -DENABLE_PERTURBATION_MODE -DENABLE_BUOYANCY_PERTURBATION_MODE"
else
   AC_MSG_RESULT([no])
fi

ENABLE_DEBUG='no'
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug], [compile in debug mode (default=no)])],
              [ENABLE_DEBUG=$enableval])

AM_CONDITIONAL([ENABLE_DEBUG], [test "$ENABLE_DEBUG" = "yes"])

AC_MSG_CHECKING([whether we are compiling in debug mode])
if test "x$ENABLE_DEBUG" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -Wall -Wno-maybe-uninitialized -Werror -g -O0"
    FCFLAGS="$FCFLAGS -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow"
else
    AC_MSG_RESULT([no])
    FCFLAGS="$FCFLAGS -O3 -funroll-all-loops -flto -DNDEBUG"
fi


ENABLE_OPENMP='no'
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--enable-openmp], [enable OpenMP (default=no)])],
              [ENABLE_OPENMP=$enableval])

AM_CONDITIONAL([ENABLE_OPENMP], [test "$ENABLE_OPENMP" = "yes"])

AC_MSG_CHECKING([whether we are enabling OpenMP])
if test "x$ENABLE_OPENMP" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -fopenmp"
else
    AC_MSG_RESULT([no])
fi


ENABLE_IW_TEST_CASE='no'
AC_ARG_ENABLE([iw-test-case],
              [AS_HELP_STRING([--enable-iw-test-case], [enable IW test case (default=no)])],
              [ENABLE_IW_TEST_CASE=$enableval])

AM_CONDITIONAL([ENABLE_IW_TEST_CASE], [test "$ENABLE_IW_TEST_CASE" = "yes"])

AC_MSG_CHECKING([whether we are enabling the IW test case])
if test "x$ENABLE_IW_TEST_CASE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_IW_TEST_CASE"
else
    AC_MSG_RESULT([no])
fi


ENABLE_RT_TEST_CASE='no'
AC_ARG_ENABLE([rt-test-case],
              [AS_HELP_STRING([--enable-rt-test-case], [enable RT test case (default=no)])],
              [ENABLE_RT_TEST_CASE=$enableval])

AM_CONDITIONAL([ENABLE_RT_TEST_CASE], [test "$ENABLE_RT_TEST_CASE" = "yes"])

AC_MSG_CHECKING([whether we are enabling the RT test case])
if test "x$ENABLE_RT_TEST_CASE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_RT_TEST_CASE"
else
    AC_MSG_RESULT([no])
fi

AC_OUTPUT
