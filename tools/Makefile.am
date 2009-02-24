AM_CPPFLAGS = -I${HEADAS}/include -I@top_srcdir@/include
AM_CFLAGS   = -Wall
#AM_CFLAGS   = -Wall -Werror

bin_PROGRAMS = binary_stream generate_photons

binary_stream_SOURCES = binary_stream.c
binary_stream_LDFLAGS  = @libcfitsio@ -lnsl
binary_stream_LDFLAGS += @libape@ -lreadline
binary_stream_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
binary_stream_LDADD = @top_srcdir@/libraries/libdetectors.a -lm

generate_photons_SOURCES  = generate_photons.c generate_photons.h
generate_photons_LDFLAGS  = @libcfitsio@ -lnsl
generate_photons_LDFLAGS += @libape@ -lreadline
generate_photons_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
generate_photons_LDADD    = @top_srcdir@/libraries/libdetectors.a
generate_photons_LDADD   += @top_srcdir@/libraries/libastrosources.a
generate_photons_LDADD   += @top_srcdir@/libraries/libsatellites.a
generate_photons_LDADD   += -lgsl -lgslcblas -lm
