AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall -Werror


########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ladspec.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ladspec

ladspec_SOURCES = ladspec.c ladspec.h
ladspec_LDADD = @top_builddir@/libraries/libsixt/libsixt.la
