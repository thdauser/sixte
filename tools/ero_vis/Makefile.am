AM_CFLAGS  = -I@top_srcdir@/simput/libsimput 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/cfitsio 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/wcslib/C
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall
#AM_CFLAGS += -Wall -Werror

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_vis.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_vis

ero_vis_SOURCES= ero_vis.c
ero_vis_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
ero_vis_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
ero_vis_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
ero_vis_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
ero_vis_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

