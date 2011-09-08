AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = runsixt.par


############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS = runsixt

runsixt_SOURCES = runsixt.c runsixt.h
runsixt_LDADD   = @top_builddir@/libraries/libsixt/libsixt.la

