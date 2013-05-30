
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
dist_pfiles_DATA=htrssim.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=htrssim

htrssim_SOURCES=htrssim.c htrssim.h
htrssim_LDADD =@top_builddir@/libsixt/libsixt.a
htrssim_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
htrssim_LDADD+=@top_builddir@/extlib/heainit/libhdinit.a
htrssim_LDADD+=@top_builddir@/extlib/heaio/libhdio.a
htrssim_LDADD+=@top_builddir@/extlib/heautils/libhdutils.a
htrssim_LDADD+=@top_builddir@/extlib/heasp/libhdsp.a
htrssim_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
htrssim_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
htrssim_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
htrssim_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a

