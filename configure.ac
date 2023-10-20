#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([sixte], [m4_esyscmd([build-aux/git-version-gen .tarball-version])], [sixte-support@lists.fau.de])

# Store auxiliary files in the sub-directory 'build-aux'.
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR(build-aux/m4)

AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_SRCDIR([libsixt/sixt.h])

# For cmocka testing with TAP output
AC_REQUIRE_AUX_FILE([tap-driver.sh])

AC_REQUIRE_AUX_FILE([git-version-gen])

# Check for programs:
AC_PROG_CC # Search for a C-Compiler and store it in the Variable 'CC'.
AC_PROG_CXX # Search for a CXX-Compiler
AC_PROG_FC # Search for a Fortran compiler.
AM_PROG_AR

AC_PROG_LIBTOOL

# find SIMPUT library
AC_ARG_WITH(simput,
	[AC_HELP_STRING([--with-simput=<SIMPUT path>],
	[Provide an alternative path to the SIMPUT software])],
	simput_dir=$withval,
	simput_dir=${prefix}
)

AC_SUBST(simput_dir)


# Compiler Flags:
# "-g -W -Wall -O0" for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# "-openmp" to use openMP multicore functionality
# "-pthread" for threading
# "-std=c++11" to support pthread
CFLAGS="$CFLAGS -g -ansi -std=c99 -pedantic -D_GNU_SOURCE -W -Wall -DHTRS_ARCPIXELS -Dg77Fortran -pthread"
CPPFLAGS="$CPPCFLAGS -I$simput_dir/include"
CXXFLAGS="$CXXFLAGS -pthread -std=c++11"

# Linker Flags:
LDFLAGS="$LDFLAGS -g -W -Wall  -L$simput_dir/lib"


OSX_VERSION_MAJOR=""
OSX_VERSION_MINOR=""
case $build_os in
  darwin*)
  OSX_VERSION_MAJOR=`sw_vers -productVersion | awk -F '.' '{print $1}'`
  OSX_VERSION_MINOR=`sw_vers -productVersion | awk -F '.' '{print $2}'`
  CFLAGS="$CFLAGS -headerpad_max_install_names"
  CPPFLAGS="$CPPFLAGS -headerpad_max_install_names"
  LDFLAGS="$LDFLAGS -headerpad_max_install_names"
  ;;
esac
AC_SUBST(OSX_VERSION_MAJOR)
AC_SUBST(OSX_VERSION_MINOR)


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

AC_SEARCH_LIBS([clock_gettime], [rt], [], [ echo "Warning:  No Real Time library found! (This is ok on a Mac OSX)"], [])

####################################
# expat
AC_SEARCH_LIBS(XML_Parse, [expat],
	[], [AC_MSG_ERROR([No expat library found!])], [])

####################################
# ncurses and readline (required by ape and progressbar library)
AC_SEARCH_LIBS(beep, ncurses, [], [AC_MSG_ERROR([No ncurses library found!])], [])
AC_SEARCH_LIBS(readline, readline, [], [AC_MSG_ERROR([No readline library found!])], [])


####################################
# GSL Libraries:
AC_ARG_WITH(gsl,
	[AC_HELP_STRING([--with-gsl=<GSL path>],
	[Deprecated (use --with-gsl-prefix instead)])],
	[AC_MSG_ERROR([Option --with-gsl deprecated. Use --with-gsl-prefix instead.])],
	[]
)


AX_PATH_GSL(1.16,[],[AC_MSG_ERROR(could not find required version of GSL >= 1.16)])

# Correct library linking information and compiler flags are stored in GSL_LIBS
# and GSL_CFLAGS after successful execution of AX_PATH_GSL.
LDFLAGS="$LDFLAGS $GSL_LIBS"
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"

AC_SEARCH_LIBS([cblas_sdot], [gslcblas], [], \
			 [AC_MSG_ERROR([No gslcblas library found!])], [])
AC_SEARCH_LIBS([gsl_sf_erf_Q], [gsl], [], \
			 [AC_MSG_ERROR([No gsl library found!])], [])


####################################
# boostlib
AX_BOOST_BASE([1.40],[],[AC_MSG_ERROR([Install with: sudo apt-get install libboost-all-dev])])

# Correct library linking information and compiler flags are stored in BOOST_LDFLAGS
# and BOOST_CPPFLAGS after successful execution of AX_BOOST_BASE.
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS"


####################################

# Check for header files:
AC_HEADER_STDC
AC_CHECK_HEADERS([limit.h math.h stdio.h stdlib.h stdarg.h string.h simput.h fitsio.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_SEARCH_LIBS(ffexist, [cfitsio], [], [AC_MSG_ERROR([ cfitsio not found (should be part of simput)!])], -lm)
AC_SEARCH_LIBS(wcssub, [wcs], [], [AC_MSG_ERROR([ libwcs not found (should be part of simput)!])], -lm)
AC_SEARCH_LIBS([fftw_free], [fftw3], [], [AC_MSG_ERROR([ libfftw not found (should be part of simput)!])], -lm)
AC_SEARCH_LIBS([ape_test], [ape], [], [AC_MSG_ERROR([ libape not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([atSun], [atFunctions], [], [AC_MSG_ERROR([ libatFunctions not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([headas_chat], [hdio], [], [AC_MSG_ERROR([ libhdio not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([posstring], [posstring], [], [AC_MSG_ERROR([ libposstring not found (should be part of simput)!])], [-lm])
AC_SEARCH_LIBS([nhinit], [labnh], [], [AC_MSG_ERROR([ liblabnh not found (should be part of simput)!])], [-lm -lposstring])
AC_SEARCH_LIBS([headas_clobberfile], [hdutils], [], [AC_MSG_ERROR([ libhdutils not found (should be part of simput)!])], [-lm -lhdio])
AC_SEARCH_LIBS([headas_init], [hdinit], [], [AC_MSG_ERROR([ libhdinit not found (should be part of simput)!])], [-lm -lhdio -lhdutils])
AC_SEARCH_LIBS([freeRMF], [simput], [], [AC_MSG_ERROR([ libsimput not found! specify path with --with-libsimput])], [-lm -lape -lhdinit -lhdio -lhdutils -lhdinit] )


# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([sixteconfig.h])

# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!

AC_CONFIG_FILES([Makefile
		libsixt/Makefile
		tools/Makefile
    tools/sixte_arfgen/Makefile
		tools/athenawfisim/Makefile
		tools/attgen_dither/Makefile
		tools/comabackpro/Makefile
		tools/comadet/Makefile
		tools/comaexp/Makefile
		tools/comaimg/Makefile
		tools/comaimgPM/Makefile
		tools/comaphovign/Makefile
		tools/comarecon/Makefile
		tools/epicmos1_events/Makefile
		tools/epicmos2_events/Makefile
		tools/epicpn_events/Makefile
		tools/ero_calevents/Makefile
		tools/ero_exposure/Makefile
		tools/ero_fits2tm/Makefile
		tools/ero_rawevents/Makefile
		tools/ero_vis/Makefile
		tools/erosim/Makefile
		tools/evpat/Makefile
		tools/exposure_map/Makefile
		tools/fudgexp/Makefile
		tools/gendetsim/Makefile
		tools/gennoisespec/Makefile
		tools/gradeddetection/Makefile
		tools/htrssim/Makefile
		tools/imgev/Makefile
		tools/ladsim/Makefile
		tools/makelc/Makefile
		tools/makespec/Makefile
		tools/nustarsim/Makefile
		tools/orbatt/Makefile
		tools/pha2pi/Makefile
		tools/phogen/Makefile
		tools/phoimg/Makefile
		tools/pixdetillum/Makefile
		tools/piximpacts/Makefile
		tools/projev/Makefile
		tools/psfgen/Makefile
		tools/pulsetemplgen/Makefile
		tools/pulsetemplimport/Makefile
		tools/radec2xy/Makefile
		tools/runsixt/Makefile
		tools/runtes/Makefile
		tools/sixteversion/Makefile
		tools/streamtotriggers/Makefile
		tools/tes_grades/Makefile
		tools/tesconstpileup/Makefile
		tools/tesreconstruction/Makefile
		tools/teslib/Makefile
		tools/tesrecons/Makefile
		tools/tessim/Makefile
		tools/tesstream/Makefile
		tools/xifupipeline/Makefile
		tools/xml2svg/Makefile
		tools/xms_pixtemp/Makefile
		tools/tesgenimpacts/Makefile
		extlib/Makefile
		extlib/progressbar/Makefile
		tools/runmask/Makefile
		test/Makefile
		test/unit/Makefile
		])

AC_OUTPUT
