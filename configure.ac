#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([sixt], [0.01], [christian.schmid@sternwarte.uni-erlangen.de])
AC_CONFIG_SRCDIR([libraries/libdetectors/wfidetector.c])
# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Check for programs:
AC_PROG_CC     # Search for a C-Compiler and store it in the Variable 'CC'.

AC_DISABLE_SHARED
LT_INIT #([disable-shared]) #([shared]) 
AC_PROG_LIBTOOL #deprecated
AM_PROG_LIBTOOL
#AC_DISABLE_SHARED
#AC_DISABLE_STATIC

# Compiler and Linker Flags:
# "-g -W -Wall -O0" for the compiler for debugging
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
CFLAGS="$CFLAGS -I$HEADAS/include -g -W -Wall -O0 -pg"
LDFLAGS="$LDFLAGS -L$HEADAS/lib -g -W -Wall -O0 -pg "


# Check for libraries:  (objdump -R/-T libncurses...)
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])
# -lnsl
#LIBS+="-lnsl"
AC_SEARCH_LIBS(ffopen, [cfitsio_3.12 cfitsio_3.07 cfitsio_3.04 cfitsio], [], \
			 [AC_MSG_ERROR([No cfitsio library found!])], [])

AC_SEARCH_LIBS(beep, ncurses, [], [AC_MSG_ERROR([No ncurses library found!])], [])
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
AC_CHECK_HEADERS([limit.h math.h malloc.h stdio.h stdlib.h stdarg.h string.h])
#AC_CHECK_HEADERS([headas.h headas_error.h]) # -> problem: not accepted by preprocessor (but by compiler)!


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST


# Checks for library functions.


# PIL parameter files have to be installed to the $PFILES directory.
AC_ARG_VAR([PFILESDIR], [target directory for .par files])
if test -z "$PFILESDIR"
then
	if test -z "$HEADAS"
	then
		PFILESDIR="\${HEADAS}/pfiles"
	else
		AC_MSG_ERROR([Environment variable '\$HEADAS' is not set], [])
	fi
fi
AC_SUBST([PFILESDIR])


# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!
AC_CONFIG_FILES([Makefile
		libraries/Makefile
		libraries/libsixt/Makefile
		libraries/libdetectors/Makefile
		libraries/libmirrors/Makefile
		libraries/libsatellites/Makefile
		libraries/libastrosources/Makefile
		tools/Makefile
		tools/photon_generation/Makefile
		tools/photon_imaging/Makefile
		tools/wfi_simulation/Makefile
		tools/framestore_simulation/Makefile
		tools/htrs_simulation/Makefile
		tools/xms_simulation/Makefile
		tools/event_projection/Makefile
		tools/eroexposure/Makefile
		tools/psfgen/Makefile
		tools/binary_stream/Makefile
		tools/htrs_psfgen/Makefile
		tools/analyse_xms_events/Makefile
		tools/pattern_recombination/Makefile])
AC_OUTPUT

