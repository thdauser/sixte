#CFLAGS = -g -W -Wall -O0       # activate this line for debug mode -> now in configure.ac

# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 

# LDADD = ...   add additional libraries


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libsixt.la

# Sources:
libsixt_la_SOURCES = sixt_random.c vector.c


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = sixt.h sixt_random.h point.h vector.h

