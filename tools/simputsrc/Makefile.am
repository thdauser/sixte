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
dist_pfiles_DATA=simputsrc.par

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=simputsrc

simputsrc_SOURCES=simputsrc.c simputsrc.h
simputsrc_LDADD =@top_builddir@/libsixt/libsixt.a
simputsrc_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.a
simputsrc_LDADD+=@top_builddir@/extlib/heainit/libhdinit.a
simputsrc_LDADD+=@top_builddir@/extlib/heaio/libhdio.a
simputsrc_LDADD+=@top_builddir@/extlib/heautils/libhdutils.a
simputsrc_LDADD+=@top_builddir@/extlib/heasp/libhdsp.a
simputsrc_LDADD+=@top_builddir@/extlib/simput/extlib/ape/src/libape.a
simputsrc_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
simputsrc_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
simputsrc_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
