AM_CPPFLAGS = -I${HEADAS}/include -I@top_srcdir@/include
#AM_CFLAGS   = -Wall -Werror
AM_CFLAGS   = -Wall

bin_PROGRAMS = generate_photons

generate_photons_SOURCES  = generate_photons.c generate_photons.h
generate_photons_LDFLAGS  = @libcfitsio@ -lnsl
generate_photons_LDFLAGS += @libape@ -lreadline
generate_photons_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
generate_photons_LDADD    = @top_srcdir@/libraries/libsatellites/libsatellites.a
generate_photons_LDADD   += @top_srcdir@/libraries/libastrosources.a
generate_photons_LDADD   += -lgsl -lgslcblas -lm
