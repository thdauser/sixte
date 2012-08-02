AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_fits2tm.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_fits2tm

ero_fits2tm_SOURCES= ero_fits2tm.c
ero_fits2tm_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
ero_fits2tm_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
ero_fits2tm_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
ero_fits2tm_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
ero_fits2tm_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la


