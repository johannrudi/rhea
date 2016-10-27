dnl
dnl This file is part of rhea2.
dnl

AC_INIT([rhea],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [johann@ices.utexas.edu])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/rhea.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/rhea_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

RHEA_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                [DEBUG])

echo "o---------------------------------------"
echo "| Checking for MPI and related programs"
echo "o---------------------------------------"

SC_MPI_CONFIG([RHEA],[yes],[yes])
SC_MPI_ENGAGE([RHEA])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
SC_C_VERSION
LT_INIT

echo "o---------------------------------------"
echo "| Checking for libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([RHEA])
P4EST_CHECK_LIBRARIES([RHEA])
MANGLL_CHECK_LIBRARIES([RHEA])
YMIR_CHECK_LIBRARIES([RHEA])
RHEA_CHECK_LIBRARIES([RHEA])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

SC_REQUIRE_FUNCS([exp log pow sqrt floor ceil erf])

echo "o---------------------------------------"
echo "| Checking third-party software"
echo "o---------------------------------------"

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

# Prepare make dist behavior and dependencies.
RHEA_DIST_DENY=
RHEA_DISTCLEAN=

SC_AS_SUBPACKAGE([RHEA])
dnl
export P4EST_SC_SOURCE=../sc
P4EST_AS_SUBPACKAGE([RHEA])
AC_DEFINE([WITH_P4EST], 1, [Use the p4est library])
dnl
export MANGLL_SC_SOURCE=../sc
export MANGLL_P4EST_SOURCE=../p4est
MANGLL_AS_SUBPACKAGE([RHEA])
dnl
export YMIR_SC_SOURCE=../sc
export YMIR_P4EST_SOURCE=../p4est
export YMIR_MANGLL_SOURCE=../mangll
YMIR_AS_SUBPACKAGE([RHEA])

# Implement make dist behavior and dependencies.
AM_CONDITIONAL([RHEA_DIST_DENY], [test "x$RHEA_DIST_DENY" != x])
AC_SUBST([RHEA_DISTCLEAN])

# Print summary.

AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
dnl AC_DEFINE_UNQUOTED(C_VERSION,   ["${C_VERSION}"],   [C compiler version])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [F77 compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(BLAS_LIBS,   ["${BLAS_LIBS}"],   [BLAS libraries])
AC_DEFINE_UNQUOTED(LAPACK_LIBS, ["${LAPACK_LIBS}"], [LAPACK libraries])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])
AC_DEFINE_UNQUOTED(FLIBS,       ["${FLIBS}"],       [Fortran libraries])

echo "o----------------------------------"
echo "| Results for rhea are"
echo "o----------------------------------"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
dnl echo "| C_VERSION:   $C_VERSION"
echo "| CFLAGS:      $CFLAGS"
echo "| F77:         $F77"
echo "| FFLAGS:      $FFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| BLAS_LIBS:   $BLAS_LIBS"
echo "| LAPACK_LIBS: $LAPACK_LIBS"
echo "| LIBS:        $LIBS"
echo "| FLIBS:       $FLIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.rhea.pre])

AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([RHEA])
P4EST_FINAL_MESSAGES([RHEA])
MANGLL_FINAL_MESSAGES([RHEA])
YMIR_FINAL_MESSAGES([RHEA])
RHEA_FINAL_MESSAGES([RHEA])
