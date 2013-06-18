AM_CFLAGS =-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
noinst_LIBRARIES=libhdsp.a

# Sources:
libhdsp_a_SOURCES=WriteRMF.c WritePHA.c WriteARF.c SPutils.c		\
                  ReturnChannel.c MergeARFRMF.c TransposeRMF.c SPio.c	\
                  RMF.c ReadRMF.c ReadARF.c ReadPHA.c

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=heaspio.h heasp.h
