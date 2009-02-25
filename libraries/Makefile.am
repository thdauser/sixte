#CFLAGS = 

AM_CPPFLAGS = -I${HEADAS}/include -I@top_srcdir@/include
AM_CFLAGS  = -Wall -Werror

lib_LIBRARIES = libastrosources.a libdetectors.a libsatellites.a

libastrosources_a_SOURCES = spectrum.c sources.c photon.c fits_pha.c random.c check_fov.c

libdetectors_a_SOURCES = detector.c eventlist.c random.c
#libdetectors_a_LDFLAGS  = @libcfitsio@ -lnsl
#libdetectors_a_LDFLAGS += @libape@ -lreadline
#libdetectors_a_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@

libsatellites_a_SOURCES = orbatt.c fits_attitude.c fits_ctlg.c vector.c

