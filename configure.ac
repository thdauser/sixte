#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([sixt], [0.01], [christian.schmid@sternwarte.uni-erlangen.de])
AC_CONFIG_SRCDIR([tools/binary_stream/binary_stream.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])



# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
#AC_DISABLE_SHARED
#AC_DISABLE_STATIC


# Checks for libraries.
CFLAGS="$CFLAGS -I$HEADAS/include"
LDFLAGS="$LDFLAGS -L$HEADAS/lib"


AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])


libcfitsio="no"
AC_CHECK_LIB(cfitsio, ffopen, [libcfitsio="-lcfitsio"], [], [-lm])
AC_CHECK_LIB(cfitsio_3.04, ffopen, [libcfitsio="-lcfitsio_3.04"], [], [-lm])
AC_CHECK_LIB(cfitsio_3.07, ffopen, [libcfitsio="-lcfitsio_3.07"], [], [-lm])
AC_CHECK_LIB(cfitsio_3.12, ffopen, [libcfitsio="-lcfitsio_3.12"], [], [-lm])
if test "x$libcfitsio" = "xno"
then
	AC_MSG_ERROR([No cfitsio library found!])
fi
AC_SUBST(libcfitsio)


libape="no"
AC_CHECK_LIB(ape_2.4.0, pilclose_, [libape="-lape_2.4.0"], [], [-lreadline])
if test "x$libape" = "xno"
then
	AC_MSG_ERROR([No ape (PIL) library found!])
fi
AC_SUBST(libape)


#libhdinit="no"
#AC_CHECK_LIB(hdinit_2.5, headas_init, [libhdinit="-lhdinit_2.5"], [], 
#  [-lhdutils_2.5 -lhdio_2.5 -lape_2.4.0 -lreadline -lcfitsio_3.12 -lnsl -lm])
#if test "x$libhdinit" = "xno"
#then
#	AC_MSG_ERROR([No hdinit library found!])
#fi
libhdinit="-lhdinit_2.5"
AC_SUBST(libhdinit)


#libhdutils="no"
#AC_CHECK_LIB(hdutils_2.5, set_toolname, [libhdutils="-lhdutils_2.5"], [], [])
#if test "x$libhdutils" = "xno"
#then
#	AC_MSG_ERROR([No hdutils library found!])
#fi
libhdutils="-lhdutils_2.5"
AC_SUBST(libhdutils)


#libhdio="no"
#AC_CHECK_LIB(hdio_2.5, headas_chat, [libhdio="-lhdio_2.5"], [], [])
#if test "x$libhdio" = "xno"
#then
#	AC_MSG_ERROR([No hdio library found!])
#fi
libhdio="-lhdio_2.5"
AC_SUBST(libhdio)



# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h])



# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST



# Checks for library functions.

AC_CONFIG_FILES([Makefile
		tools/Makefile
		tools/binary_stream/Makefile
		libraries/Makefile])
#		libraries/libdetectors/Makefile])
AC_OUTPUT
