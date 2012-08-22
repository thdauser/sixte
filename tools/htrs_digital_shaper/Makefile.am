AM_CFLAGS =-I@top_srcdir@/libraries/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-Wall 
#AM_CFLAGS+=-Wall -Werror

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=htrs_digital_shaper.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=htrs_digital_shaper

htrs_digital_shaper_SOURCES=htrs_digital_shaper.c htrs_digital_shaper.h
htrs_digital_shaper_LDADD =@top_builddir@/libraries/libsixt/libsixt.la
htrs_digital_shaper_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.la 
htrs_digital_shaper_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
htrs_digital_shaper_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/libfftw3.la


