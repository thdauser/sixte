AM_CFLAGS =-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libhdinit.la

# Sources:
libhdinit_la_SOURCES=headas_init.c

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=headas.h
