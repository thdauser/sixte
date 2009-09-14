AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 
#AM_CFLAGS += -DNORMALIZE_RMF # explicitly normalize the RSP to an RMF


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = erosita.eventlist.tpl 
dist_fitstemplates_DATA	+= wfi.full1024.eventlist.tpl 
dist_fitstemplates_DATA += wfi.window16.eventlist.tpl
dist_fitstemplates_DATA += xms.eventlist.tpl
dist_fitstemplates_DATA += htrs.eventlist.tpl


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libdetectors.la

# Sources:
libdetectors_la_SOURCES = genericdetector.c squarepixels.c hexagonalpixels.c \
			  eventfile.c impactlist.c \
			  simplecosmicbackground.c \
			  framestoredetector.c erositaeventfile.c \
			  wfidetector.c wfieventfile.c \
			  xmsdetector.c xmseventfile.c \
			  htrsdetector.c htrseventfile.c

libdetectors_la_LIBADD  = @top_builddir@/libraries/libsixt/libsixt.la
libdetectors_la_LIBADD += @top_builddir@/libraries/libastrosources/libastrosources.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = genericdetector.h squarepixels.h hexagonalpixels.h \
		  eventfile.h impactlist.h \
		  simplecosmicbackground.h \
		  framestoredetector.h erositaeventfile.h erositaevent.h \
		  wfidetector.h wfieventfile.h wfievent.h \
		  xmsdetector.h xmseventfile.h xmsevent.h \
		  htrsdetector.h htrseventfile.h htrsevent.h

