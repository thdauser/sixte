#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([sixt], [0.01], [christian.schmid@sternwarte.uni-erlangen.de])
AC_CONFIG_SRCDIR([src/binary_stream.c])
# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for programs.
AC_PROG_CC     # Search for a C-Compiler and store it in the Variable 'CC'.
AC_DISABLE_SHARED
#AC_DISABLE_STATIC
AC_PROG_LIBTOOL


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


#AC_CHECK_LIB(readline, ?, [], AC_MSG_ERROR([No readline library found!]), [])
#AC_CHECK_LIB(ncurses, ?, [], AC_MSG_ERROR([No ncurses library found!]), [])

libape="no"
AC_CHECK_LIB(ape_2.4.0, pilclose_, [libape="-lape_2.4.0"], [], [-lreadline -lncurses])
if test "x$libape" = "xno"
then
	AC_MSG_ERROR([No ape (PIL) library found!])
fi
AC_SUBST(libape)


libhdinit="no"
#AC_CHECK_LIB(hdinit_2.5, headas_init, [libhdinit="-lhdinit_2.5"], [], 
#  [-lhdutils_2.5 -lhdio_2.5 -lape_2.4.0 -lreadline -lncurses -lcfitsio_3.12 -lnsl -lm])
#if test "x$libhdinit" = "xno"
#then
#	AC_MSG_ERROR([No hdinit library found!])
#fi
libhdinit="-lhdinit_2.5"
AC_SUBST(libhdinit)


libhdutils="no"
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



# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!
AC_CONFIG_FILES([Makefile
		src/Makefile])
AC_OUTPUT
