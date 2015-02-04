AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=epicpn_events.par

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=epicpn_events

epicpn_events_SOURCES=epicpn_events.c epicpn_events.h
epicpn_events_LDADD =@top_builddir@/libsixt/libsixt.la
