AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = analyse_xms_events.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = analyse_xms_events

analyse_xms_events_SOURCES= analyse_xms_events.c analyse_xms_events.h
analyse_xms_events_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
analyse_xms_events_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
analyse_xms_events_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
analyse_xms_events_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
analyse_xms_events_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

