dnl
dnl This file is part of ymir.
dnl

AC_INIT([ymir],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [ccgo@ices.utexas.edu])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/ymir.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/ymir_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

YMIR_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                [DEBUG])
YMIR_ARG_DISABLE([vtk-doubles], [use doubles for vtk file data],
                 [VTK_DOUBLES])
YMIR_ARG_DISABLE([vtk-binary], [write vtk ascii file data],
                 [VTK_BINARY])
YMIR_ARG_DISABLE([vtk-zlib], [disable zlib compression for vtk binary data],
                 [VTK_COMPRESSION])

YMIR_ARG_DISABLE([component-checks],
                 [disable checks for subcomponents (libsc, p4est, etc.)],
                 [COMPONENT_CHECKS])
YMIR_ARG_DISABLE([examples], [disable examples],[EXAMPLES])

echo "o---------------------------------------"
echo "| Checking for MPI and related programs"
echo "o---------------------------------------"

SC_MPI_CONFIG([YMIR],[yes],[yes])
SC_MPI_ENGAGE([YMIR])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
SC_C_VERSION
LT_INIT

echo "o---------------------------------------"
echo "| Checking for libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([YMIR])
P4EST_CHECK_LIBRARIES([YMIR])
MANGLL_CHECK_LIBRARIES([YMIR])
YMIR_CHECK_LIBRARIES([YMIR])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

SC_REQUIRE_FUNCS([memset pow sqrt floor fsync])
SC_REQUIRE_FUNCS([wordexp wordfree])

echo "o---------------------------------------"
echo "| Checking third-party software"
echo "o---------------------------------------"

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

# Prepare make dist behavior and dependencies.
YMIR_DIST_DENY=
YMIR_DISTCLEAN=

SC_AS_SUBPACKAGE([YMIR])
export P4EST_SC_SOURCE=../sc
P4EST_AS_SUBPACKAGE([YMIR])
AC_DEFINE([WITH_P4EST], 1, [Use the p4est library])
export MANGLL_SC_SOURCE=../sc
export MANGLL_P4EST_SOURCE=../p4est
MANGLL_AS_SUBPACKAGE([YMIR])
MANGLL_ENABLE_EXAMPLES="no" # disable mangll examples TODO: not working!

# Implement make dist behavior and dependencies.
AM_CONDITIONAL([YMIR_DIST_DENY], [test "x$YMIR_DIST_DENY" != x])
AC_SUBST([YMIR_DISTCLEAN])

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
echo "| Results for ymir are"
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
AC_CONFIG_FILES([Makefile Makefile.ymir.pre])
YMIR_CONFIG_SCRIPT([example/performance/ymir_mass_vec_convergence_plot.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_trilinear_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_mass_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_derivative_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_pressure_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_vel_rob_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_velocity_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_stokes_elem.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_submesh.sh])
YMIR_CONFIG_SCRIPT([test/ymir_test_mass_vec.sh])

AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([YMIR])
P4EST_FINAL_MESSAGES([YMIR])
MANGLL_FINAL_MESSAGES([YMIR])
YMIR_FINAL_MESSAGES([YMIR])
