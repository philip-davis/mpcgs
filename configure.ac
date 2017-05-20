AC_INIT([mpcgs], [0.2.0])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

AC_PROG_CC([nvcc])
AC_PROG_CXX([nvcc])

AC_DEFINE(SFMT_MEXP,19937)

LIBS+="-lm"

AC_OUTPUT(Makefile)
