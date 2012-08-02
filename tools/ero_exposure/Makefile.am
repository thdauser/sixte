AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall
#AM_CFLAGS += -Wall -Werror


########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_exposure.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_exposure

ero_exposure_SOURCES = ero_exposure.c
ero_exposure_LDADD = @top_builddir@/libraries/libsixt/libsixt.la
ero_exposure_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
ero_exposure_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
ero_exposure_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
ero_exposure_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

