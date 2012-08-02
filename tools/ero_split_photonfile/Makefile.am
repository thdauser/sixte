AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_split_photonfile.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_split_photonfile

ero_split_photonfile_SOURCES = ero_split_photonfile.c ero_split_photonfile.h
ero_split_photonfile_LDADD   = @top_builddir@/libraries/libsixt/libsixt.la
ero_split_photonfile_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
ero_split_photonfile_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
ero_split_photonfile_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
ero_split_photonfile_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

