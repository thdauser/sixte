# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = erosita_pattern_recombination.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = erosita_pattern_recombination

erosita_pattern_recombination_SOURCES = erosita_pattern_recombination.c erosita_pattern_recombination.h
erosita_pattern_recombination_LDADD  = @top_builddir@/libraries/libdetectors/libdetectors.la

