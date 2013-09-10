AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heasp
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=htrs_pixel_intensity.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=htrs_pixel_intensity

htrs_pixel_intensity_SOURCES=htrs_pixel_intensity.c htrs_pixel_intensity.h
htrs_pixel_intensity_LDADD =@top_builddir@/libsixt/libsixt.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/heainit/libhdinit.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/heaio/libhdio.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/heautils/libhdutils.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/heasp/libhdsp.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
htrs_pixel_intensity_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
