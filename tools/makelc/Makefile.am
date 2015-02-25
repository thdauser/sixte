AM_CPPFLAGS =-I@top_srcdir@/libsixt
 

########## DIRECTORIES ##############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=makelc.par

AM_CFLAGS=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=makelc

makelc_SOURCES=makelc.c makelc.h
makelc_LDADD =@top_builddir@/libsixt/libsixt.la
