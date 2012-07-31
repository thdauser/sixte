AM_CFLAGS = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS+= -Wall
#AM_CFLAGS += -Wall -Werror

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = ero_vis.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = ero_vis

ero_vis_SOURCES = ero_vis.c
ero_vis_LDADD = @top_builddir@/libraries/libsixt/libsixt.la

