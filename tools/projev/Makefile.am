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

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=projev.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=projev

projev_SOURCES=projev.c
projev_LDADD =@top_builddir@/libsixt/libsixt.a
projev_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/heainit/libhdinit.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/heaio/libhdio.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/heautils/libhdutils.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/heasp/libhdsp.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
projev_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
