AM_CFLAGS  = -I@top_srcdir@/simput/libsimput 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/cfitsio 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/wcslib/C
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = runsixt.par


############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS = runsixt

runsixt_SOURCES= runsixt.c runsixt.h
runsixt_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
runsixt_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
runsixt_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
runsixt_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
runsixt_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

