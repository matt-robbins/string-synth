AC_INIT([String-Synth], 1.0)
AM_INIT_AUTOMAKE
AC_PROG_CC

AC_HEADER_STDC
AC_CHECK_HEADERS([jack/jack.h math.h])

AC_CONFIG_SRCDIR([src/string-synth.c])
AC_CONFIG_HEADER([src/config.h])

AC_CHECK_LIB([m], [sin], [], [AC_MSG_ERROR([string-synth requires libm])])
AC_CHECK_LIB([jack], [jack_client_open], [], \
    [AC_MSG_ERROR([string-synth requires JACK])])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
