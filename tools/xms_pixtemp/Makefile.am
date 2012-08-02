AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = xms_pixtemp.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = xms_pixtemp

xms_pixtemp_SOURCES = xms_pixtemp.c xms_pixtemp.h
xms_pixtemp_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
xms_pixtemp_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
xms_pixtemp_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
xms_pixtemp_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
xms_pixtemp_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

