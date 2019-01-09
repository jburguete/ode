AC_PREREQ([2.69])
AC_INIT([ode], [0.0.0], [https://github.com/jburguete/ode])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC(mpicc, egcc, gcc, clang)
AC_MSG_NOTICE([Compiler: $CC])
AC_LANG([C])
AC_PROG_MAKE_SET
PKG_PROG_PKG_CONFIG

# Checks for libraries.
AC_CHECK_LIB([m], [sqrtl])
PKG_CHECK_MODULES([GSL], [gsl])
PKG_CHECK_MODULES([XML], [libxml-2.0])
PKG_CHECK_MODULES([GLIB], [glib-2.0])
PKG_CHECK_MODULES([GTHREAD], [gthread-2.0])

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([float.h libintl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([sqrtl])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
