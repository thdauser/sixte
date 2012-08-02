AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = comaimg.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = comaimg

comaimg_SOURCES= comaimg.c comaimg.h
comaimg_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
comaimg_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
comaimg_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
comaimg_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
comaimg_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la


