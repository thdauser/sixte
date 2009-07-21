AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/src # REMOVE THIS LINE !!
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libdetectors.la

# Sources:
libdetectors_la_SOURCES = detector.c \
			  squarepixels.c genericdetector.c \
			  eventfile.c impactlist.c \
			  detectorbackground.c \
			  framestoredetector.c erositaeventfile.c \
			  wfidetector.c wfieventfile.c
libdetectors_la_LIBADD 	= @top_srcdir@/libraries/libsixt/libsixt.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = detector.h detectors.types.h detectors.def.h detectors.enum.h \
		  genericdetector.h squarepixels.h \
		  eventfile.h impactlist.h \
		  detectorbackground.h \
		  framestoredetector.h erositaeventfile.h erositaevent.h \
		  wfidetector.h wfieventfile.h wfievent.h

