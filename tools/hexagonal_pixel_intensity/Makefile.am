#CFLAGS = -g -W -Wall -O0       # activate this line for debug mode -> now in configure.ac

# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
AM_CFLAGS += -I@top_srcdir@/libraries/libmirrors
AM_CFLAGS += -I@top_srcdir@/libraries/libsatellites
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
#dist_pfiles_DATA = hexagonal_pixel_intensity.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = hexagonal_pixel_intensity

hexagonal_pixel_intensity_SOURCES = hexagonal_pixel_intensity.c
hexagonal_pixel_intensity_LDADD  = @top_builddir@/libraries/libdetectors/libdetectors.la
hexagonal_pixel_intensity_LDADD += @top_builddir@/libraries/libmirrors/libmirrors.la

