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

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=testbkg

testbkg_SOURCES=testbkg.c
testbkg_LDADD =@top_builddir@/libsixt/libsixt.a
testbkg_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a 
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/heainit/libhdinit.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/heaio/libhdio.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/heautils/libhdutils.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/heasp/libhdsp.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
testbkg_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
