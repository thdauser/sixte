AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/progressbar/include


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=sixte_arfgen.par

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=sixte_arfgen

sixte_arfgen_SOURCES=sixte_arfgen.c sixte_arfgen.h
sixte_arfgen_LDADD =@top_builddir@/libsixt/libsixt.la
