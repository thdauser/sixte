CFLAGS = -g -W -Wall -O0       # activate this line for debug mode

# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include -I@top_srcdir@/include
AM_CFLAGS += -Wall
#AM_CFLAGS += -Wall -Werror

# LDADD = ...   add additional libraries


pfilesdir   = @PFILESDIR@
pfiles_DATA = 

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES = libastrosources.a libdetectors.a libsatellites.a

# Sources:
libastrosources_a_SOURCES = check_fov.c fits_pha.c pointsources.c photon.c \
			random.c spectrum.c
libdetectors_a_SOURCES = detector.c eventlist.c psf.c vignetting.c random.c vector.c
libsatellites_a_SOURCES = fits_attitude.c fits_ctlg.c orbatt.c vector.c


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = binary_stream htrs_psfgen psfgen \
		photon_generation photon_detection photon_imaging event_projection

binary_stream_SOURCES  = binary_stream.c
binary_stream_LDADD    = libdetectors.a
pfiles_DATA           += binary_stream.par

htrs_psfgen_SOURCES  = htrs_psfgen.c
htrs_psfgen_LDADD    = libdetectors.a

photon_generation_SOURCES  = photon_generation.c
photon_generation_LDADD    = libdetectors.a
photon_generation_LDADD   += libastrosources.a
photon_generation_LDADD   += libsatellites.a
pfiles_DATA               += photon_generation.par

photon_detection_SOURCES  = photon_detection.c
photon_detection_LDADD    = libdetectors.a
photon_detection_LDADD   += libastrosources.a
photon_detection_LDADD   += libsatellites.a
pfiles_DATA              += photon_detection.par

photon_imaging_SOURCES  = photon_imaging.c
photon_imaging_LDADD    = libdetectors.a
photon_imaging_LDADD   += libastrosources.a
photon_imaging_LDADD   += libsatellites.a
pfiles_DATA            += photon_imaging.par

event_projection_SOURCES  = event_projection.c
event_projection_LDADD    = libdetectors.a
event_projection_LDADD   += libastrosources.a
event_projection_LDADD   += libsatellites.a
pfiles_DATA              += event_projection.par

psfgen_SOURCES  = psfgen.c
psfgen_LDADD    = libdetectors.a
psfgen_LDADD   += libastrosources.a
psfgen_LDADD   += libsatellites.a
pfiles_DATA    += psfgen.par

