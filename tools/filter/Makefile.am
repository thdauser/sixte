AM_CXXFLAGS =-I@top_srcdir@/libsixt
AM_CXXFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio
AM_CXXFLAGS+=-I@top_srcdir@/extlib/Utils
AM_CXXFLAGS+=-I@top_srcdir@/extlib/simput/extlib/fftw/api
AM_CXXFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
#pfilesdir=$(pkgdatadir)/pfiles
#dist_pfiles_DATA=filter.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=filter

filter_SOURCES=filter.cpp filter.h
filter_LDADD =@top_builddir@/libsixt/libsixt.a
filter_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
filter_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
filter_LDADD+=@top_builddir@/extlib/Utils/libutilsic.a
