AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 

# Use a constant light curve instead of red noise according to Timmer & Koenig (1995):
AM_CFLAGS += -DCONSTANT_LIGHTCURVE


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libastrosources.la

# Sources:
libastrosources_la_SOURCES = pointsources.c sourceimage.c photon.c spectrum.c \
			fits_pha.c check_fov.c
libastrosources_la_LIBADD = @top_builddir@/libraries/libsixt/libsixt.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = astrosources.h pointsources.h sourceimage.h photon.h spectrum.h \
		fits_pha.h check_fov.h
