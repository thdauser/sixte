AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/progressbar/include 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=ero_fits2tm.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=ero_fits2tm

ero_fits2tm_SOURCES=ero_fits2tm.c
ero_fits2tm_LDADD =@top_builddir@/libsixt/libsixt.la
