AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = wfmphovign.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = wfmphovign

wfmphovign_SOURCES= wfmphovign.c wfmphovign.h
wfmphovign_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
wfmphovign_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
wfmphovign_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
wfmphovign_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
wfmphovign_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la
