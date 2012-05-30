AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_gti.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_gti

ero_gti_SOURCES = ero_gti.c
ero_gti_LDADD = @top_builddir@/libraries/libsixt/libsixt.la

