AM_CFLAGS  = -I${HEADAS}/include -I@top_srcdir@/include
AM_CFLAGS += -Wall
#AM_CFLAGS += -Wall -Werror
CFLAGS = -g -W -Wall -O0       # activate this line for debug mode

#AM_LDFLAGS = -L${HEADAS}/lib
LDADD  = @libhdinit@ @libhdutils@ @libhdio@
LDADD += @libape@ -lreadline -lncurses 
LDADD += @libcfitsio@ -lnsl -lm


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES = libastrosources.a libdetectors.a libsatellites.a

# Sources:
libastrosources_a_SOURCES = spectrum.c sources.c pointsources.c photon.c fits_pha.c random.c check_fov.c
libdetectors_a_SOURCES = detector.c eventlist.c psf.c random.c
libsatellites_a_SOURCES = orbatt.c fits_attitude.c fits_ctlg.c vector.c


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = binary_stream photon_generation photon_detection photon_imaging psfgen

binary_stream_SOURCES  = binary_stream.c
binary_stream_LDADD    = ${LDADD} libdetectors.a

photon_generation_SOURCES  = photon_generation.c
photon_generation_LDADD    = ${LDADD} -lgsl -lgslcblas
photon_generation_LDADD   += libdetectors.a
photon_generation_LDADD   += libastrosources.a
photon_generation_LDADD   += libsatellites.a

photon_detection_SOURCES  = photon_detection.c
photon_detection_LDADD    = ${LDADD} -lgsl -lgslcblas
photon_detection_LDADD   += libdetectors.a
photon_detection_LDADD   += libastrosources.a
photon_detection_LDADD   += libsatellites.a

photon_imaging_SOURCES  = photon_imaging.c
photon_imaging_LDADD    = ${LDADD} -lgsl -lgslcblas
photon_imaging_LDADD   += libdetectors.a
photon_imaging_LDADD   += libastrosources.a
photon_imaging_LDADD   += libsatellites.a

psfgen_SOURCES  = psfgen.c
psfgen_LDADD    = ${LDADD} -lgsl -lgslcblas
psfgen_LDADD   += libdetectors.a
psfgen_LDADD   += libastrosources.a
psfgen_LDADD   += libsatellites.a

