#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([sixt], [0.15], [christian.schmid@sternwarte.uni-erlangen.de])
AC_CONFIG_AUX_DIR([build-aux]) # Store auxiliary files in the sub-directory 'build-aux'.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([libraries/libsixt/sixt.h])

# Check if HEADAS is available.
if test -z "$HEADAS"; then
   AC_MSG_WARN([HEADAS is not defined - the build might fail!])
fi

# Check for programs:
AC_PROG_CC # Search for a C-Compiler and store it in the Variable 'CC'.

LT_INIT #([disable-shared]) #([shared]) 

# Compiler Flags:
# "-g -W -Wall -O0" for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# "-openmp" to use openMP multicore functionality
CFLAGS="$CFLAGS -I$HEADAS/include -g -fopenmp -ansi -std=c99 -pedantic -D_GNU_SOURCE -W -Wall -DHTRS_ARCPIXELS -DCORRECT_RMF_FIRST_CHANNEL"

# Linker Flags:
LDFLAGS="$LDFLAGS -L$HEADAS/lib -g -W -Wall"

# Check for libraries:  (objdump -R/-T libncurses...)
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])

####################################
# cfitsio
AC_SEARCH_LIBS(ffopen, 
	[cfitsio_3.27 cfitsio_3.25 cfitsio_3.24 cfitsio_3.22 cfitsio_3.18 cfitsio_3.12 cfitsio_3.07 cfitsio_3.04 cfitsio], 
	[], [AC_MSG_ERROR([No supported version of cfitsio library found!])], [])

####################################
# expat
AC_SEARCH_LIBS(XML_Parse, [expat], 
	[], [AC_MSG_ERROR([No expat library found!])], [])

#################################### 
# ncurses and readline (required by ape library)
AC_SEARCH_LIBS(beep, ncurses, [], [AC_MSG_ERROR([No ncurses library found!])], [])
AC_SEARCH_LIBS(readline, readline, [], [AC_MSG_ERROR([No readline library found!])], [])

####################################
# PIL / APE
lib_ape="no"
AC_MSG_CHECKING([for versioned libape in \$HEADAS])
my_ape=`ls "${HEADAS}"/lib/libape*.so 2> /dev/null`
if ! test -z "$my_ape"
then
   my_ape=`basename "$my_ape"`
   my_ape=`expr "$my_ape" : 'lib\(ape_.*\)\.so'`
fi
if ! test -z "$my_ape"
then
   AC_MSG_RESULT([yes ($my_ape)])
else
   AC_MSG_RESULT([no])
fi
AC_SEARCH_LIBS(ape_trad_init, ape $my_ape, [],
          [AC_MSG_ERROR([No ape library found!])], [])
AC_SUBST(lib_ape)

# HEADAS libraries:
####################################
# HDIO
lib_hdio="no"
AC_MSG_CHECKING([for versioned hdio in \$HEADAS])
my_hdio=`ls "${HEADAS}"/lib/libhdio*.so 2> /dev/null`
if ! test -z "$my_hdio"
then
   my_hdio=`basename "$my_hdio"`
   my_hdio=`expr "$my_hdio" : 'lib\(hdio_.*\)\.so'`
fi
if ! test -z "$my_hdio"
then
   AC_MSG_RESULT([yes ($my_hdio)])
else
   AC_MSG_RESULT([no])
fi
AC_SEARCH_LIBS(headas_chat, hdio $my_hdio, [],
          [AC_MSG_ERROR([No hdio library found!])], [])
AC_SUBST(lib_hdio)
####################################
# HDUTILS
lib_hdutils="no"
AC_MSG_CHECKING([for versioned hdutils in \$HEADAS])
my_hdutils=`ls "${HEADAS}"/lib/libhdutils*.so 2> /dev/null`
if ! test -z "$my_hdutils"
then
   my_hdutils=`basename "$my_hdutils"`
   my_hdutils=`expr "$my_hdutils" : 'lib\(hdutils_.*\)\.so'`
fi
if ! test -z "$my_hdutils"
then
   AC_MSG_RESULT([yes ($my_hdutils)])
else
   AC_MSG_RESULT([no])
fi
AC_SEARCH_LIBS(set_toolname, hdutils $my_hdutils, [],
          [AC_MSG_ERROR([No hdutils library found!])], [])
AC_SUBST(lib_hdutils)
#################################### 
# HDINIT
lib_hdinit="no"
AC_MSG_CHECKING([for versioned hdinit in \$HEADAS])
my_hdinit=`ls "${HEADAS}"/lib/libhdinit*.so 2> /dev/null`
if ! test -z "$my_hdinit"
then
   my_hdinit=`basename "$my_hdinit"`
   my_hdinit=`expr "$my_hdinit" : 'lib\(hdinit_.*\)\.so'`
fi
if ! test -z "$my_hdinit"
then
   AC_MSG_RESULT([yes ($my_hdinit)])
else
   AC_MSG_RESULT([no])
fi
AC_SEARCH_LIBS(headas_init, hdinit $my_hdinit, [],
          [AC_MSG_ERROR([No hdinit library found!])], [])
AC_SUBST(lib_hdinit)
#################################### 
# HDSP
lib_hdsp="no"
AC_MSG_CHECKING([for versioned hdsp in \$HEADAS])
my_hdsp=`ls "${HEADAS}"/lib/libhdsp*.so 2> /dev/null`
if ! test -z "$my_hdsp"
then
   my_hdsp=`basename "$my_hdsp"`
   my_hdsp=`expr "$my_hdsp" : 'lib\(hdsp_.*\)\.so'`
fi
if ! test -z "$my_hdsp"
then
   AC_MSG_RESULT([yes ($my_hdsp)])
else
   AC_MSG_RESULT([no])
fi
AC_SEARCH_LIBS(ReadRMFMatrix, hdsp $my_hdsp, [],
              [AC_MSG_ERROR([No hdsp library found!])], [])
AC_SUBST(lib_hdsp)
#################################### 
# WCS Library:
AC_SEARCH_LIBS(wcsbth, [wcs_6.11 wcs_6.10 wcs_6.9 wcs-4.4.4 wcs-4.3.3 wcs-4.3], [], \
		       [AC_MSG_ERROR([No WCS library found!])], [])

#################################### 
# GSL Libraries:
AC_SEARCH_LIBS(cblas_sdot, [gslcblas], [], \
			 [AC_MSG_ERROR([No gslcblas library found!])], [])
AC_SEARCH_LIBS(gsl_sf_erf_Q, [gsl], [], \
			 [AC_MSG_ERROR([No gsl library found!])], [])


####################################
# SIMPUT
AC_SEARCH_LIBS(getSimputSource, [simput], [], 
	       [AC_MSG_ERROR([No supported version of SIMPUT library found!])], [])

#################################### 
# Check for header files:
AC_HEADER_STDC
AC_CHECK_HEADERS([limit.h math.h malloc.h stdio.h stdlib.h stdarg.h string.h])
#AC_CHECK_HEADERS([headas.h headas_error.h]) 
# -> problem: not accepted by preprocessor (but by compiler)!


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
		libraries/Makefile
		libraries/libsixt/Makefile
		tools/Makefile
		tools/runsixt/Makefile
		tools/phogen/Makefile
		tools/phoimg/Makefile
		tools/phovign/Makefile
		tools/ladsim/Makefile
		tools/ladspec/Makefile
		tools/genlc/Makefile
		tools/psfgen/Makefile
		tools/xms_pixtemp/Makefile
		tools/ero_exposure/Makefile
		tools/ero_events/Makefile
		tools/ero_fits2tm/Makefile
		tools/comaexp/Makefile
		tools/ero_split_photonfile/Makefile
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
		tools/mmap/Makefile
		share/Makefile
		share/instruments/Makefile])
AC_OUTPUT

