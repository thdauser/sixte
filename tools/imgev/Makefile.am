AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-Wall 

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=imgev.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=imgev

imgev_SOURCES=imgev.c
imgev_LDADD =@top_builddir@/libsixt/libsixt.a
imgev_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
imgev_LDADD+=@top_builddir@/extlib/heainit/libhdinit.a
imgev_LDADD+=@top_builddir@/extlib/heaio/libhdio.a
imgev_LDADD+=@top_builddir@/extlib/heautils/libhdutils.a
imgev_LDADD+=@top_builddir@/extlib/heasp/libhdsp.a
imgev_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
imgev_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
imgev_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
imgev_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
