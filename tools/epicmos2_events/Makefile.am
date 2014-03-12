AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/heasp
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/fftw/api
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=epicmos2_events.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=epicmos2_events

epicmos2_events_SOURCES=epicmos2_events.c epicmos2_events.h
epicmos2_events_LDADD =@top_builddir@/libsixt/libsixt.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/heainit/libhdinit.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/heaio/libhdio.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/heautils/libhdutils.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/heasp/libhdsp.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
epicmos2_events_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
