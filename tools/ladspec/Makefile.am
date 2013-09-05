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

########## DIRECTORIES ##############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=ladspec.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=ladspec

ladspec_SOURCES=ladspec.c ladspec.h
ladspec_LDADD =@top_builddir@/libsixt/libsixt.a
ladspec_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/heainit/libhdinit.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/heaio/libhdio.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/heautils/libhdutils.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/heasp/libhdsp.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
ladspec_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
