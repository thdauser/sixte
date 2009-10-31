# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libmirrors
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -I@top_srcdir@/libraries/libsatellites
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = comarecon.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = comarecon

comarecon_SOURCES = comarecon.c comarecon.h
comarecon_LDADD  = @top_builddir@/libraries/libmirrors/libmirrors.la
comarecon_LDADD += @top_builddir@/libraries/libastrosources/libastrosources.la
comarecon_LDADD += @top_builddir@/libraries/libsatellites/libsatellites.la
comarecon_LDADD += @top_builddir@/libraries/libdetectors/libdetectors.la

