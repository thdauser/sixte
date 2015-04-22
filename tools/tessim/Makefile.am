AM_CPPFLAGS =-I@top_srcdir@/libsixt
 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=tessim.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=tessim

tessim_SOURCES=tessim.c tessim.h tessim_datastream.c tessim_impactlist.c tessim_tesrecord.c tes_simulation.c
tessim_LDADD =@top_builddir@/libsixt/libsixt.la
