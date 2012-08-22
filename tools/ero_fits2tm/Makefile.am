AM_CFLAGS =-I@top_srcdir@/libraries/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-Wall 
#AM_CFLAGS+=-Wall -Werror

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=ero_fits2tm.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=ero_fits2tm

ero_fits2tm_SOURCES=ero_fits2tm.c
ero_fits2tm_LDADD =@top_builddir@/libraries/libsixt/libsixt.la
ero_fits2tm_LDADD+=@top_builddir@/extlib/simput/libsimput/libsimput.la 
ero_fits2tm_LDADD+=@top_builddir@/extlib/simput/extlib/wcslib/C/libwcs-4.13.4.a
ero_fits2tm_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/libfftw3.la
