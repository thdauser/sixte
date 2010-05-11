# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_split_photonfile.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_split_photonfile

ero_split_photonfile_SOURCES = ero_split_photonfile.c ero_split_photonfile.h
ero_split_photonfile_LDADD   = @top_builddir@/libraries/libastrosources/libastrosources.la

