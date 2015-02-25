AM_CPPFLAGS =-I@top_srcdir@/libsixt
 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=xms_pixtemp.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=xms_pixtemp

xms_pixtemp_SOURCES=xms_pixtemp.c xms_pixtemp.h
xms_pixtemp_LDADD =@top_builddir@/libsixt/libsixt.la

