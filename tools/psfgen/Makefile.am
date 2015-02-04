AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=psfgen.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=psfgen

psfgen_SOURCES=psfgen.c
psfgen_LDADD =@top_builddir@/libsixt/libsixt.la
