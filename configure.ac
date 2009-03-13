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


# Check for programs:
AC_PROG_CC     # Search for a C-Compiler and store it in the Variable 'CC'.

AC_DISABLE_SHARED
#AC_DISABLE_STATIC
AC_PROG_LIBTOOL

# Compiler and Linker Flags:
CFLAGS="$CFLAGS -I$HEADAS/include"
LDFLAGS="$LDFLAGS -L$HEADAS/lib"


# Check for libraries:
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])
# -lnsl
AC_SEARCH_LIBS(ffopen, [cfitsio_3.12 cfitsio_3.07 cfitsio_3.04 cfitsio], [], \
			 [AC_MSG_ERROR([No cfitsio library found!])], [])

# -lncurses (beep ?)     (objdump -R/-T libncurses... )
#LIBS+="-lncurses -lnsl"
AC_SEARCH_LIBS(readline, readline, [], [AC_MSG_ERROR([No readline library found!])], [])
AC_SEARCH_LIBS(pilclose_, ape_2.4.0, [], \
			 [AC_MSG_ERROR([No ape library (PIL) found!])], [])

# HEADAS libraries:
AC_SEARCH_LIBS(headas_chat, [hdio_2.5], [], \
			 [AC_MSG_ERROR([No hdio library found!])], [])
AC_SEARCH_LIBS(set_toolname, [hdutils_2.5], [], \
			 [AC_MSG_ERROR([No hdutils library found!])], [])
AC_SEARCH_LIBS(headas_init, [hdinit_2.5], [], \
			 [AC_MSG_ERROR([No hdinit library found!])], [])
AC_SEARCH_LIBS(ReadRMFMatrix, [hdsp_2.5], [], \
			 [AC_MSG_ERROR([No hdsp library found!])], [])

# GSL Libraries:
AC_SEARCH_LIBS(cblas_sdot, [gslcblas], [], \
			 [AC_MSG_ERROR([No gslcblas library found!])], [])
AC_SEARCH_LIBS(gsl_sf_erf_Q, [gsl], [], \
			 [AC_MSG_ERROR([No gsl library found!])], [])


# Check for header files:
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

