AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/progressbar/include

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=teslib.par

############ BINARIES ################# 

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=teslib

teslib_SOURCES=teslib.c teslib.h
teslib_LDADD =@top_builddir@/libsixt/libsixt.la
