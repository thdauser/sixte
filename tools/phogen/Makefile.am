AM_CFLAGS =-I@includedir@
AM_CFLAGS +=-I@top_srcdir@/libsixt
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=phogen.par

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=phogen

phogen_SOURCES=phogen.c phogen.h
phogen_LDADD =@top_builddir@/libsixt/libsixt.la
