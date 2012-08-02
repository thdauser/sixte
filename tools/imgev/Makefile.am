AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 


########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = imgev.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = imgev

imgev_SOURCES= imgev.c
imgev_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
imgev_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
imgev_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
imgev_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
imgev_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la
