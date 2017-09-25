AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/progressbar/include 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=ero_calevents.par

AM_CFLAGS=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=ero_calevents

ero_calevents_SOURCES=ero_calevents.c ero_calevents.h
ero_calevents_LDADD =@top_builddir@/libsixt/libsixt.la
