AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/Utils
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
#pfilesdir=$(pkgdatadir)/pfiles
#dist_pfiles_DATA=trigger.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=trigger

trigger_SOURCES=trigger.cpp trigger.h
trigger_LDADD =@top_builddir@/libsixt/libsixt.a
trigger_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/heainit/libhdinit.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/heaio/libhdio.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/heautils/libhdutils.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/heasp/libhdsp.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
trigger_LDADD+=@top_builddir@/extlib/simput/extlib/Utils/libutilsic.a
