#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([sixte], [2.0.1], [thomas.dauser@sternwarte.uni-erlangen.de])

# Store auxiliary files in the sub-directory 'build-aux'.
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([libsixt/sixt.h])

# Check for programs:
AC_PROG_CC # Search for a C-Compiler and store it in the Variable 'CC'.
AC_PROG_CXX # Search for a CXX-Compiler 
AC_PROG_FC # Search for a Fortran compiler.
AM_PROG_AR

AC_PROG_LIBTOOL

AC_LIB_HAVE_LINKFLAGS([simput])

# Compiler Flags:
# "-g -W -Wall -O0" for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# "-openmp" to use openMP multicore functionality
CFLAGS="$CFLAGS -g -ansi -std=c99 -pedantic -D_GNU_SOURCE -W -Wall -DHTRS_ARCPIXELS -Dg77Fortran"

# Linker Flags:
# LDFLAGS="$LDFLAGS -g -W -Wall"
LDFLAGS="$LDFLAGS -g -W -Wall -L${prefix}/lib -L$LIBSIMPUT_PREFIX/lib"

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
# boostlib
AC_CONFIG_MACRO_DIR(./build-aux/m4/)
AX_BOOST_BASE([1.40],[],[AC_MSG_ERROR([Install with: sudo apt-get install libboost-all-dev])])

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
AC_SEARCH_LIBS([freeRMF], [simput], [], [AC_MSG_ERROR([ libsimput not found! specify path with --with-libsimput-prefix])], [-lm -lape -lhdinit -lhdio -lhdutils -llabnh -lposstring -lhdinit] )
AC_SEARCH_LIBS(ffexist, [cfitsio], [], [AC_MSG_ERROR([ cfitsio not found (should be part of simput)!])], -lm)
AC_SEARCH_LIBS(wcssub, [wcs], [], [AC_MSG_ERROR([ libwcs not found (should be part of simput)!])], -lm)
AC_SEARCH_LIBS([fftw_free], [fftw3], [], [AC_MSG_ERROR([ libfftw not found (should be part of simput)!])], -lm)
AC_SEARCH_LIBS([ape_test], [ape], [], [AC_MSG_ERROR([ libape not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([atSun], [atFunctions], [], [AC_MSG_ERROR([ libatFunctions not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([headas_chat], [hdio], [], [AC_MSG_ERROR([ libhdio not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([headas_clobberfile], [hdutils], [], [AC_MSG_ERROR([ libhdutils not found (should be part of simput)!])], [-lm -lhdio])
AC_SEARCH_LIBS([headas_init], [hdinit], [], [AC_MSG_ERROR([ libhdinit not found (should be part of simput)!])], [-lm -lhdio -lhdutils])
AC_SEARCH_LIBS([posstring], [posstring], [], [AC_MSG_ERROR([ libposstring not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([nhinit], [labnh], [], [AC_MSG_ERROR([ liblabnh not found (should be part of simput)!])], [-lm -lposstring])


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
		tools/makespec/Makefile
		tools/makelc/Makefile
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
		tools/htrssim/Makefile
		tools/tes_grades/Makefile
		tools/gendetsim/Makefile
		tools/comaimg/Makefile
		tools/comaimgPM/Makefile
		tools/comadet/Makefile
		tools/comarecon/Makefile
		tools/comabackpro/Makefile
		tools/projev/Makefile
		tools/evpat/Makefile
		tools/imgev/Makefile
		tools/fudgexp/Makefile
		tools/orbatt/Makefile
		tools/epicpn_events/Makefile
		tools/epicmos1_events/Makefile
		tools/epicmos2_events/Makefile
		tools/piximpacts/Makefile
		tools/tesstream/Makefile
		tools/pixdetillum/Makefile
		tools/pulsetemplgen/Makefile
		tools/pulsetemplimport/Makefile
		tools/streamtotriggers/Makefile
		tools/runtes/Makefile
		tools/trigger/Makefile
		tools/pulsegrade/Makefile
		tools/filter/Makefile
		tools/pseudoenergy/Makefile
		tools/energy/Makefile
		tools/gennoisespec/Makefile
		tools/tesconstpileup/Makefile
		extlib/Makefile
		extlib/Utils/Makefile])

AC_OUTPUT

