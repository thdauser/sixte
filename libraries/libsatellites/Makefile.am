AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/libraries/libsatellites
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libsatellites.la

# Sources:
libsatellites_la_LIBADD  = @top_builddir@/libraries/libsixt/libsixt.la
libsatellites_la_SOURCES = attitudecatalog.c attitudefile.c


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = attitudecatalog.h attitudefile.h telescope.h
