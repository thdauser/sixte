# If CFLAGS is overwritten, the -I${HEADAS}/include flag has to be reset manually:
AM_CFLAGS  = -I${HEADAS}/include 
AM_CFLAGS += -I@top_srcdir@/libraries/libsixt 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

#LDADD = ...    add additional libraries


########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir = @PFILESDIR@
dist_pfiles_DATA = comaimg.par


############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS = comaimg

comaimg_SOURCES = comaimg.c comaimg.h
comaimg_LDADD  = @top_builddir@/libraries/libsixt/libsixt.la


