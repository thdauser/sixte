AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = comarecon.par

AM_CFLAGS += -DSIXT_DATA_PATH='"$(pkgdatadir)"'


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = comarecon

comarecon_SOURCES= comarecon.c comarecon.h
comarecon_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
comarecon_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
comarecon_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
comarecon_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
comarecon_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la



