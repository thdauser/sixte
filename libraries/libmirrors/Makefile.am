AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -I@top_srcdir@/libraries/libsatellites
AM_CFLAGS += -I@top_srcdir@/libraries/libmirrors
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 
#AM_CFLAGS += -DNORMALIZE_RMF # explicitly normalize the RSP to an RMF


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libmirrors.la

# Sources:
libmirrors_la_SOURCES = psf.c vignetting.c codedmask.c
libmirrors_la_LIBADD  = @top_builddir@/libraries/libsixt/libsixt.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = psf.h vignetting.h codedmask.h
