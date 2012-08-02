AM_CFLAGS  = -I@top_srcdir@/simput/libsimput 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/cfitsio 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/wcslib/C
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ladsim.par


############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS = ladsim

ladsim_SOURCES= ladsim.c ladsim.h
ladsim_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
ladsim_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
ladsim_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
ladsim_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
ladsim_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la
