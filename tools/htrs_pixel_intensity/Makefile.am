AM_CFLAGS  = -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = htrs_pixel_intensity.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = htrs_pixel_intensity

htrs_pixel_intensity_SOURCES= htrs_pixel_intensity.c htrs_pixel_intensity.h
htrs_pixel_intensity_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la
htrs_pixel_intensity_LDADD += @top_builddir@/simput/libsimput/libsimput.la 
htrs_pixel_intensity_LDADD += @top_builddir@/simput/extlib/cfitsio/libcfitsio.a 
htrs_pixel_intensity_LDADD += @top_builddir@/simput/extlib/wcslib/C/libwcs-4.13.4.a
htrs_pixel_intensity_LDADD += @top_builddir@/simput/extlib/fftw/libfftw3.la

