AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=analyse_xms_events.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=analyse_xms_events

analyse_xms_events_SOURCES=analyse_xms_events.c analyse_xms_events.h
analyse_xms_events_LDADD =@top_builddir@/libsixt/libsixt.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/heainit/libhdinit.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/heaio/libhdio.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/heautils/libhdutils.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/heasp/libhdsp.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
analyse_xms_events_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a

