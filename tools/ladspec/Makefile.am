AM_CFLAGS  = -I@top_srcdir@/simput/libsimput 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/cfitsio 
AM_CFLAGS += -I@top_srcdir@/simput/extlib/wcslib/C
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ladspec.par

AM_CFLAGS += -DSIXT_DATA_PATH='"$(pkgdatadir)"'


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ladspec

ladspec_SOURCES= ladspec.c ladspec.h
ladspec_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
ladspec_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
ladspec_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
ladspec_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
ladspec_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la
