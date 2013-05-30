AM_CFLAGS =-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
noinst_LIBRARIES=libhdinit.a

# Sources:
libhdinit_a_SOURCES=headas_init.c

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=headas.h
