AM_CPPFLAGS =-I@top_srcdir@/libsixt
 

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=exposure_map.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=exposure_map

exposure_map_SOURCES=exposure_map.c
exposure_map_LDADD =@top_builddir@/libsixt/libsixt.la
