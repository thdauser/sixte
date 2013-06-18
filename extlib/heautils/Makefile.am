AM_CFLAGS =-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/ape/include

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
noinst_LIBRARIES=libhdutils.a

# Sources:
libhdutils_a_SOURCES=headas_sort.c headas_lock.c headas_history.c	\
                     headas_file_check.c headas_error_init.c		\
                     HDgtcalf_test.c headas_utils.c headas_toolname.c	\
                     headas_svdfit.c headas_smooth.c headas_rand.c	\
                     headas_gti.c headas_error_map.c headas_error.c	\
                     HDgtcalf_f77_wrap.c HDgtcalf.c mt.c		\
                     headas_polyfit.c headas_parstamp.c			\
                     headas_error_manager.c headas_copykeys.c		\
                     HDgtcalf_standalone.c HDgtcalf_internal.c

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=headas_gti.h mt.h headas_rand.h headas_polyfit.h	\
                headas_error_cfitsio.h HDgtcalf_standalone.h		\
                HDgtcalf_internal.h headas_utils.h headas_lock.h	\
                headas_error_internal.h headas_error.h			\
                headas_copykeys.h hdcal.h
