
dnl RHEA_CHECK_RHEAKIT(PREFIX)
dnl Check for the RHEAKIT library
dnl

AC_DEFUN([RHEA_CHECK_RHEAKIT], [

AC_MSG_CHECKING([for rheakit])

SC_ARG_WITH_PREFIX([rheakit], [set <dir> to rhea-kit installation], [RHEAKIT],
                   [$1], [=<dir>])
if test "x$$1_WITH_RHEAKIT" != "xno" ; then
  if test "x$$1_WITH_RHEAKIT" = "xyes" ; then
    $1_RHEAKIT_INCLUDES=""
    $1_RHEAKIT_LIBS=""
  else
    AC_MSG_CHECKING([Rhea-kit include directory and Makefiles])
    $1_RHEAKIT_DIR="$$1_WITH_RHEAKIT"
    if test ! -d "$$1_RHEAKIT_DIR" ; then
      AC_MSG_ERROR([$$1_RHEAKIT_DIR not found])
    fi
    $1_RHEAKIT_INC_DIR="$$1_RHEAKIT_DIR/include"
    if test ! -d "$$1_RHEAKIT_INC_DIR" ; then
      AC_MSG_ERROR([$$1_RHEAKIT_INC_DIR not found])
    fi
    $1_RHEAKIT_LIB_DIR="$$1_RHEAKIT_DIR/lib"
    if test ! -d "$$1_RHEAKIT_LIB_DIR" ; then
      AC_MSG_ERROR([$$1_RHEAKIT_LIB_DIR not found])
    fi
    $1_RHEAKIT_INCLUDES="-I$$1_RHEAKIT_INC_DIR"
    $1_RHEAKIT_LIBS="-L$$1_RHEAKIT_LIB_DIR -lblock"
  fi
  AC_SUBST([$1_RHEAKIT_INCLUDES])
  AC_SUBST([$1_RHEAKIT_LIBS])
  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi
])
