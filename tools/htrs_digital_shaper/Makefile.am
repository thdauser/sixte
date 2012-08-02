AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = htrs_digital_shaper.par

AM_CFLAGS += -DSIXT_DATA_PATH='"$(pkgdatadir)"'


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = htrs_digital_shaper

htrs_digital_shaper_SOURCES= htrs_digital_shaper.c htrs_digital_shaper.h
htrs_digital_shaper_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
htrs_digital_shaper_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
htrs_digital_shaper_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
htrs_digital_shaper_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
htrs_digital_shaper_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la


