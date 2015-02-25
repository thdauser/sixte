AM_CPPFLAGS =-I@top_srcdir@/libsixt
 

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=ero_exposure.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=ero_exposure

ero_exposure_SOURCES=ero_exposure.c
ero_exposure_LDADD =@top_builddir@/libsixt/libsixt.la
