#
# Configure wrapper for the progressbar library
#
# Joern Wilms, 2015.06.20
#

AUTOMAKE_OPTIONS=subdir-objects

# libraries that are installed in the lib directory
lib_LTLIBRARIES=libprogressbar.la

# include options
libprogressbar_la_CFLAGS=-I./include

# Sources
libprogressbar_la_SOURCES=lib/progressbar.c lib/statusbar.c

# Headers that wil be installed in $(includedir):
include_HEADERS=include/progressbar.h include/statusbar.h

