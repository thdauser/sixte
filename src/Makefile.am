#CFLAGS = -g -W -Wall -O0       # activate this line for debug mode -> now in configure.ac

# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include -I@top_srcdir@/include
AM_CFLAGS += -Wall 
#AM_CFLAGS += -DNORMALIZE_RMF # explicitly normalize the RSP to an RMF
AM_CFLAGS += -DCONSTANT_LIGHTCURVE # use a constant light curve instead of red noise according to Timmer & Koenig (1995)
#AM_CFLAGS += -Wall -Werror 

# LDADD = ...   add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir   = @PFILESDIR@
dist_pfiles_DATA = 

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = erosita.eventlist.tpl 
dist_fitstemplates_DATA	+= wfi.full1024.eventlist.tpl 
dist_fitstemplates_DATA += wfi.window16.eventlist.tpl


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES = libastrosources.a \
		libdetectors.a \
		libmirror.a libsatellites.a libfitsfiles.a

# Sources:
libastrosources_a_SOURCES = check_fov.c fits_pha.c pointsources.c photon.c \
			    random_sixt.c spectrum.c sourceimage.c
libdetectors_a_SOURCES    = detector.c random_sixt.c vector.c \
			squarepixels.c genericdetector.c eventfile.c detectorbackground.c \
			framestoredetector.c erositaeventfile.c \
			wfidetector.c wfieventfile.c
libmirror_a_SOURCES       = psf.c vignetting.c vector.c random_sixt.c
libsatellites_a_SOURCES   = attitudecatalog.c attitudefile.c vector.c
libfitsfiles_a_SOURCES    = impactlist.c


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = 	binary_stream htrs_psfgen psfgen \
	       	photon_generation photon_imaging event_projection \
		framestore_simulation wfi_simulation #photon_detection

binary_stream_SOURCES  = binary_stream.c
binary_stream_LDADD    = libdetectors.a libfitsfiles.a
binary_stream_LDADD   += libfitsfiles.a
dist_pfiles_DATA           += binary_stream.par

htrs_psfgen_SOURCES  = htrs_psfgen.c
htrs_psfgen_LDADD    = libmirror.a libdetectors.a libfitsfiles.a

photon_generation_SOURCES  = photon_generation.c photon_generation.h
photon_generation_LDADD    = libdetectors.a libfitsfiles.a
photon_generation_LDADD   += libastrosources.a
photon_generation_LDADD   += libsatellites.a
dist_pfiles_DATA               += photon_generation.par

#photon_detection_SOURCES  = photon_detection.c
#photon_detection_LDADD    = libdetectors.a
#photon_detection_LDADD   += libastrosources.a
#photon_detection_LDADD   += libsatellites.a
#photon_detection_LDADD   += libfitsfiles.a
#dist_pfiles_DATA              += photon_detection.par

framestore_simulation_SOURCES  = framestore_simulation.c framestore_simulation.h
framestore_simulation_LDADD    = libdetectors.a
framestore_simulation_LDADD   += libastrosources.a
framestore_simulation_LDADD   += libsatellites.a
framestore_simulation_LDADD   += libfitsfiles.a
dist_pfiles_DATA    += framestore_simulation.par

wfi_simulation_SOURCES  = wfi_simulation.c wfi_simulation.h
wfi_simulation_LDADD    = libdetectors.a
wfi_simulation_LDADD   += libastrosources.a
wfi_simulation_LDADD   += libsatellites.a
wfi_simulation_LDADD   += libfitsfiles.a
dist_pfiles_DATA    += wfi_simulation.par

photon_imaging_SOURCES  = photon_imaging.c photon_imaging.h
photon_imaging_LDADD    = libmirror.a
photon_imaging_LDADD   += libastrosources.a
photon_imaging_LDADD   += libsatellites.a
dist_pfiles_DATA            += photon_imaging.par

event_projection_SOURCES  = event_projection.c
event_projection_LDADD    = libdetectors.a
event_projection_LDADD   += libastrosources.a
event_projection_LDADD   += libsatellites.a
event_projection_LDADD   += libfitsfiles.a
dist_pfiles_DATA              += event_projection.par

psfgen_SOURCES  = psfgen.c
psfgen_LDADD    = libmirror.a
psfgen_LDADD   += libastrosources.a
psfgen_LDADD   += libsatellites.a
dist_pfiles_DATA    += psfgen.par

