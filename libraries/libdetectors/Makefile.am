AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/src # REMOVE THIS LINE !!
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 
#AM_CFLAGS += -DNORMALIZE_RMF # explicitly normalize the RSP to an RMF


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = erosita.eventlist.tpl 
dist_fitstemplates_DATA	+= wfi.full1024.eventlist.tpl 
dist_fitstemplates_DATA += wfi.window16.eventlist.tpl


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libdetectors.la

# Sources:
libdetectors_la_SOURCES = squarepixels.c genericdetector.c \
			  eventfile.c impactlist.c \
			  detectorbackground.c \
			  framestoredetector.c erositaeventfile.c \
			  wfidetector.c wfieventfile.c \
			  xmsdetector.c xmseventfile.c
libdetectors_la_LIBADD = @top_builddir@/libraries/libsixt/libsixt.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = genericdetector.h squarepixels.h \
		  eventfile.h impactlist.h \
		  detectorbackground.h \
		  framestoredetector.h erositaeventfile.h erositaevent.h \
		  wfidetector.h wfieventfile.h wfievent.h \
		  xmsdetector.h xmseventfile.h xmsevent.h

