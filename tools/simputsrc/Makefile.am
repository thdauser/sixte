AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = simputsrc.par


############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS = simputsrc

simputsrc_SOURCES = simputsrc.c simputsrc.h
simputsrc_LDADD   = @top_builddir@/libraries/libsixt/libsixt.la

