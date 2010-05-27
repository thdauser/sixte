AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt
AM_CFLAGS += -I@top_srcdir@/libraries/libastrosources
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = photonlist.tpl 


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libastrosources.la

# Sources:
libastrosources_la_SOURCES = pointsources.c extendedsources.c		\
			     sourceimage.c photon.c spectrum.c		\
			     fits_pha.c check_fov.c photonlistfile.c	\
			     linlightcurve.c kdtree.c sourcelist.c
libastrosources_la_LIBADD = @top_builddir@/libraries/libsixt/libsixt.la


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = astrosources.h pointsources.h extendedsources.h	\
		  sourceimage.h photon.h spectrum.h fits_pha.h		\
		  check_fov.h photonlistfile.h linlightcurve.h		\
		  kdtree.h sourcelist.h

