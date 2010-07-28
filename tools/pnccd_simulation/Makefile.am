# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
#AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = pnccd_simulation.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = 	pnccd_simulation

pnccd_simulation_SOURCES = pnccd_simulation.c pnccd_simulation.h
pnccd_simulation_LDADD = @top_builddir@/libraries/libdetectors/libdetectors.la

