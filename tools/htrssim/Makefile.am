AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall -Werror

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = htrssim.par

AM_CFLAGS += -DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = htrssim

htrssim_SOURCES = htrssim.c htrssim.h
htrssim_LDADD = @top_builddir@/libraries/libsixt/libsixt.la
htrssim_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
htrssim_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
htrssim_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
htrssim_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

