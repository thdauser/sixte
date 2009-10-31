#CFLAGS = -g -W -Wall -O0       # activate this line for debug mode -> now in configure.ac

# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
#AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -I@top_srcdir@/libraries/libsatellites
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = photon_generation.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = photon_generation

photon_generation_SOURCES = photon_generation.c photon_generation.h
photon_generation_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
photon_generation_LDADD += @top_builddir@/libraries/libdetectors/libdetectors.la
photon_generation_LDADD += @top_builddir@/libraries/libastrosources/libastrosources.la
photon_generation_LDADD += @top_builddir@/libraries/libsatellites/libsatellites.la

