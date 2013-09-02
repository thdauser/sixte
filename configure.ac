#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([sixte], [0.94], [christian.schmid@sternwarte.uni-erlangen.de])
AC_CONFIG_AUX_DIR([build-aux]) # Store auxiliary files in the sub-directory 'build-aux'.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([libsixt/sixt.h])

# Check for programs:
AC_PROG_CC # Search for a C-Compiler and store it in the Variable 'CC'.
AC_PROG_FC # Search for a Fortran compiler.
AC_PROG_RANLIB

# Compiler Flags:
# "-g -W -Wall -O0" for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# "-openmp" to use openMP multicore functionality
CFLAGS="$CFLAGS -g -fopenmp -ansi -std=c99 -pedantic -D_GNU_SOURCE -W -Wall -DHTRS_ARCPIXELS -Dg77Fortran"

# Linker Flags:
LDFLAGS="$LDFLAGS -g -W -Wall"

####################################
# Conditional usage of RCL library (Remeis random number server).
AC_ARG_ENABLE([rcl], 
              [AS_HELP_STRING([--enable-rcl], 
	          [enable usage of RCL library])],
	          [rcl=true], [rcl=false])
AM_CONDITIONAL([RCL], [test x$rcl = xtrue])

####################################
# Check for libraries: (objdump -R/-T libncurses...)
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])

####################################
# expat
AC_SEARCH_LIBS(XML_Parse, [expat], 
	[], [AC_MSG_ERROR([No expat library found!])], [])

#################################### 
# ncurses and readline (required by ape library)
#AC_SEARCH_LIBS(beep, ncurses, [], [AC_MSG_ERROR([No ncurses library found!])], [])
#AC_SEARCH_LIBS(readline, readline, [], [AC_MSG_ERROR([No readline library found!])], [])

#################################### 
# GSL Libraries:
AC_SEARCH_LIBS(cblas_sdot, [gslcblas], [], \
			 [AC_MSG_ERROR([No gslcblas library found!])], [])
AC_SEARCH_LIBS(gsl_sf_erf_Q, [gsl], [], \
			 [AC_MSG_ERROR([No gsl library found!])], [])

#################################### 
# Check for header files:
AC_HEADER_STDC
AC_CHECK_HEADERS([limit.h math.h stdio.h stdlib.h stdarg.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.

# PIL parameter files have to be installed to the $PFILES directory.
AC_ARG_VAR([PFILESDIR], [target directory for .par files])
if test -z "$PFILESDIR"
then
	if test -z "$HEADAS"
	then
		AC_MSG_ERROR([Environment variable '\$HEADAS' is not set], [])
	else
		PFILESDIR="\${HEADAS}/pfiles"
	fi
fi
AC_SUBST([PFILESDIR])

# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([config.h])

# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!
AC_CONFIG_FILES([Makefile
		libsixt/Makefile
		tools/Makefile
		tools/runsixt/Makefile
		tools/phogen/Makefile
		tools/phoimg/Makefile
		tools/comaphovign/Makefile
		tools/ladsim/Makefile
		tools/ladspec/Makefile
		tools/genlc/Makefile
		tools/psfgen/Makefile
		tools/xms_pixtemp/Makefile
		tools/ero_exposure/Makefile
		tools/ero_vis/Makefile
		tools/ero_calevents/Makefile
		tools/ero_rawevents/Makefile
		tools/ero_fits2tm/Makefile
		tools/erosim/Makefile
		tools/comaexp/Makefile
		tools/nustarsim/Makefile
		tools/athenapwfisim/Makefile
		tools/photonmerge/Makefile
		tools/htrssim/Makefile
		tools/htrs_psfgen/Makefile
		tools/htrs_pixel_intensity/Makefile
		tools/htrs_fits2tm/Makefile
		tools/htrs_digital_shaper/Makefile
		tools/analyse_xms_events/Makefile
		tools/gendetsim/Makefile
		tools/comaimg/Makefile
		tools/comadet/Makefile
		tools/comarecon/Makefile
		tools/projev/Makefile
		tools/evpat/Makefile
		tools/imgev/Makefile
		tools/simputmerge/Makefile
		tools/simputsrc/Makefile
		tools/fudgexp/Makefile
		tools/orbatt/Makefile
		tools/mmap/Makefile
		share/Makefile
		share/instruments/Makefile
		extlib/Makefile
		extlib/heautils/Makefile
		extlib/heaio/Makefile
		extlib/heainit/Makefile
		extlib/heasp/Makefile])

AC_CONFIG_SUBDIRS([extlib/simput])

AC_OUTPUT

