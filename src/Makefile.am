#CFLAGS = -g -W -Wall -O0       # activate this line for debug mode -> now in configure.ac

# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -I@top_srcdir@/libraries/libdetectors
AM_CFLAGS += -Wall 
#AM_CFLAGS += -DNORMALIZE_RMF # explicitly normalize the RSP to an RMF
AM_CFLAGS += -DCONSTANT_LIGHTCURVE # use a constant light curve instead of red noise according to Timmer & Koenig (1995)
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = 


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libastrosources.la \
		libmirror.la libsatellites.la

# Sources:
libastrosources_la_SOURCES = check_fov.c fits_pha.c pointsources.c photon.c \
			    spectrum.c sourceimage.c
libastrosources_la_LIBADD  = @top_builddir@/libraries/libsixt/libsixt.la

#libdetectors_la_SOURCES = detector.c \
#			  squarepixels.c genericdetector.c \
#			  eventfile.c impactlist.c \
#			  detectorbackground.c \
#			  framestoredetector.c erositaeventfile.c \
#			  wfidetector.c wfieventfile.c
#libdetectors_la_LIBADD 	= @top_builddir@/libraries/libsixt/libsixt.la

libmirror_la_SOURCES    = psf.c vignetting.c
libmirror_la_LIBADD	= @top_builddir@/libraries/libsixt/libsixt.la

libsatellites_la_SOURCES   = attitudecatalog.c attitudefile.c
libsatellites_la_LIBADD    = @top_builddir@/libraries/libsixt/libsixt.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = telescope.h check_fov.h \
		spectrum.h photon.h astrosources.h pointsources.h sourceimage.h \
		fits_pha.h attitudefile.h attitudecatalog.h psf.h vignetting.h


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = 	binary_stream htrs_psfgen psfgen #photon_detection

binary_stream_SOURCES = binary_stream.c
binary_stream_LDADD   = @top_builddir@/libraries/libdetectors/libdetectors.la
dist_pfiles_DATA     += binary_stream.par

htrs_psfgen_SOURCES  = htrs_psfgen.c
htrs_psfgen_LDADD    = libmirror.la @top_builddir@/libraries/libdetectors/libdetectors.la

#photon_detection_SOURCES  = photon_detection.c
#photon_detection_LDADD    = @top_builddir@/libraries/libdetectors/libdetectors.la
#photon_detection_LDADD   += libastrosources.la
#photon_detection_LDADD   += libsatellites.la
#dist_pfiles_DATA         += photon_detection.par

psfgen_SOURCES  = psfgen.c
psfgen_LDADD    = libmirror.la
psfgen_LDADD   += libastrosources.la
psfgen_LDADD   += libsatellites.la
dist_pfiles_DATA    += psfgen.par

