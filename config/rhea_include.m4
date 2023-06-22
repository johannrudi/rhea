dnl
dnl RHEA acinclude.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl RHEA_ARG_ENABLE(NAME, COMMENT, TOKEN, [HELPEXTRA])
dnl Check for --enable/disable-NAME using shell variable RHEA_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional RHEA_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([RHEA_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [RHEA], [$4])])

dnl RHEA_ARG_DISABLE(NAME, COMMENT, TOKEN, [HELPEXTRA])
dnl Check for --enable/disable-NAME using shell variable RHEA_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional RHEA_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([RHEA_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [RHEA], [$4])])

dnl RHEA_ARG_WITH(NAME, COMMENT, TOKEN, [HELPEXTRA])
dnl Check for --with/without-NAME using shell variable RHEA_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional RHEA_TOKEN
dnl Default is without
dnl
AC_DEFUN([RHEA_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [RHEA], [$4])])

dnl RHEA_ARG_WITHOUT(NAME, COMMENT, TOKEN, [HELPEXTRA])
dnl Check for --with/without-NAME using shell variable RHEA_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional RHEA_TOKEN
dnl Default is with
dnl
AC_DEFUN([RHEA_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [RHEA], [$4])])

dnl RHEA_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by rhea.  It can be used by other packages that
dnl link to rhea to add appropriate options to LIBS.
dnl
AC_DEFUN([RHEA_CHECK_LIBRARIES],
[
dnl SC_REQUIRE_LIB([m], [fabs])
RHEA_CHECK_RHEAKIT([$1])
])

dnl RHEA_AS_SUBPACKAGE(PREFIX, prefix)
dnl Call from a package that is using rhea as a subpackage.
dnl Sets PREFIX_DIST_DENY="yes" if rhea is make install'd.
dnl
AC_DEFUN([RHEA_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1],[m4_tolower([$1])],[RHEA],[rhea])])

dnl RHEA_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([RHEA_FINAL_MESSAGES],)

dnl RHEA_CONFIG_SCRIPT (SCRIPT)
dnl mark script as config file and set executable bit
dnl
AC_DEFUN([RHEA_CONFIG_SCRIPT],
[
AC_CONFIG_FILES([$1], [chmod +x "$1"])
])
