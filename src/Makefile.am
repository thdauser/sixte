AM_CPPFLAGS = -I${HEADAS}/include -I@top_srcdir@/include
AM_CFLAGS   = -Wall
#AM_CFLAGS   = -Wall -Werror
#CFLAGS = 


# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = binary_stream generate_photons

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES = libastrosources.a libdetectors.a libsatellites.a

# Sources:
binary_stream_SOURCES = binary_stream.c
binary_stream_LDFLAGS  = @libcfitsio@ -lnsl
binary_stream_LDFLAGS += @libape@ -lreadline
binary_stream_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
binary_stream_LDADD = @top_srcdir@/libraries/libdetectors.a -lm

generate_photons_SOURCES  = generate_photons.c generate_photons.h
generate_photons_LDFLAGS  = @libcfitsio@ -lnsl
generate_photons_LDFLAGS += @libape@ -lreadline
generate_photons_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
generate_photons_LDFLAGS += -lgsl -lgslcblas -lm
generate_photons_LDADD    = @top_srcdir@/libraries/libdetectors.a
generate_photons_LDADD   += @top_srcdir@/libraries/libastrosources.a
generate_photons_LDADD   += @top_srcdir@/libraries/libsatellites.a

libastrosources_a_SOURCES = spectrum.c sources.c photon.c fits_pha.c random.c check_fov.c

libdetectors_a_SOURCES = detector.c eventlist.c random.c
#libdetectors_a_LDFLAGS  = @libcfitsio@ -lnsl
#libdetectors_a_LDFLAGS += @libape@ -lreadline
#libdetectors_a_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@

libsatellites_a_SOURCES = orbatt.c fits_attitude.c fits_ctlg.c vector.c

