AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=epicmos1_events.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=epicmos1_events

epicmos1_events_SOURCES=epicmos1_events.c epicmos1_events.h
epicmos1_events_LDADD =@top_builddir@/libsixt/libsixt.la
