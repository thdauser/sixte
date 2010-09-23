# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -Wall #-fopenmp
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_exposure.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_exposure

ero_exposure_SOURCES = ero_exposure.c
ero_exposure_LDADD = @top_builddir@/libraries/libdetectors/libdetectors.la
ero_exposure_LDADD+= @top_builddir@/libraries/libastrosources/libastrosources.la

