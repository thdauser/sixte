AM_CPPFLAGS =-I@top_srcdir@/libsixt

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=tesgenimpacts.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=tesgenimpacts

tesgenimpacts_SOURCES=tesgenimpacts.c tesgenimpacts.h 
tesgenimpacts_LDADD =@top_builddir@/libsixt/libsixt.la