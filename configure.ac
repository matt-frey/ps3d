AC_INIT([ps3d], [0.0.1], [mf248@st-andrews.ac.uk], [], [https://github.com/matt-frey/ps3d])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_FC([gfortran])
AC_LANG(Fortran)
AM_PROG_AR
LT_INIT

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)

FCFLAGS="-std=f2008 -fdefault-real-8 -fdefault-double-8 -cpp -fPIC -mcmodel=large"

AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([src/utils/config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    src/utils/Makefile
    src/netcdf/Makefile
    src/fft/Makefile
    unit-tests/Makefile
])

#######################################################################################
##
## "--with" flags
##

## define some default search paths
declare default_lib_paths="/lib:/usr/lib:/lib64:/usr/lib64:/usr/local/lib:$LD_LIBRARY_PATH"
declare default_include_paths="/usr/include:/usr/local/lib:$C_INCLUDE_PATH:$CPLUS_INCLUDE_PATH"

# 11 March 2021
# https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=<path>], [path to netcdf library @<:@default=$NETCDF_DIR@:>@])],
            [NETCDF_DIR=$withval],
            [])

if test -n "${NETCDF_DIR}"; then
    FCFLAGS="$FCFLAGS -I$NETCDF_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_DIR/lib"
    LIBS="$LIBS -lnetcdff -lnetcdf"
else
    for path in ${default_lib_paths//:/ }; do
        AC_MSG_CHECKING([for netcdf Fortran libraries in $path])
        if test -e "$path/libnetcdff.a"; then
            NETCDF_LIBRARY_DIR="$path"
            AC_MSG_RESULT([yes])
            break
        fi
        AC_MSG_RESULT([no])
    done

    for path in ${default_include_paths//:/ }; do
        AC_MSG_CHECKING([for netcdf Fortran modules in $path])
        if test -e "$path/netcdf.mod"; then
            NETCDF_INCLUDE_DIR="$path"
            AC_MSG_RESULT([yes])
            break
        fi
        AC_MSG_RESULT([no])
    done

    if test -n "${NETCDF_LIBRARY_DIR}" && test -n "${NETCDF_INCLUDE_DIR}"; then
        FCFLAGS="$FCFLAGS -I$NETCDF_INCLUDE_DIR"
        LDFLAGS="$LDFLAGS -L$NETCDF_LIBRARY_DIR"
        LIBS="$LIBS -lnetcdff -lnetcdf"
    else
        AC_MSG_ERROR([Cannot find netcdf.
                    Please export the environment variable NETCDF_DIR pointing
                    to the root directory of netcdf.])
    fi
fi

# 11 March 2021
# https://github.com/PhysicsofFluids/AFiD/blob/master/configure.ac#L139
AC_MSG_CHECKING([whether we can compile a netcdf program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([], [
    use netcdf
    implicit none])],
[netcdf_found=yes],
[netcdf_found=no AC_MSG_ERROR([Cannot compile a netcdf program])])
AC_MSG_RESULT([$netcdf_found])

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

AC_MSG_CHECKING([whether we are compiling in verbose mode])
if test "x$ENABLE_BUOYANCY" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_BUOYANCY"
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

AC_OUTPUT
