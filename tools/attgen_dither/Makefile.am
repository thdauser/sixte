AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/progressbar/include 

########## PIL .par-FILES##############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=attgen_dither.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=attgen_dither

attgen_dither_SOURCES=attgen_dither.c
attgen_dither_LDADD =@top_builddir@/libsixt/libsixt.la
