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
dist_pfiles_DATA=athenapwfisim.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=athenapwfisim

athenapwfisim_SOURCES=athenapwfisim.c athenapwfisim.h
athenapwfisim_LDADD =@top_builddir@/libsixt/libsixt.a
athenapwfisim_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
athenapwfisim_LDADD+=@top_builddir@/extlib/heainit/libhdinit.a
athenapwfisim_LDADD+=@top_builddir@/extlib/heaio/libhdio.a
athenapwfisim_LDADD+=@top_builddir@/extlib/heautils/libhdutils.a
athenapwfisim_LDADD+=@top_builddir@/extlib/heasp/libhdsp.a
athenapwfisim_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
athenapwfisim_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
athenapwfisim_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
athenapwfisim_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
