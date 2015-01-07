AM_CXXFLAGS =-I@top_srcdir@/libsixt
AM_CXXFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio
AM_CXXFLAGS+=-I@top_srcdir@/extlib/Utils
AM_CXXFLAGS+=-I@top_srcdir@/extlib/simput/extlib/fftw/api
AM_CXXFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
#pfilesdir=$(pkgdatadir)/pfiles
#dist_pfiles_DATA=energy.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=energy

energy_SOURCES=energy.cpp energy.h
energy_LDADD =@top_builddir@/libsixt/libsixt.a
energy_LDADD+=@top_builddir@/extlib/simput/extlib/fftw/.libs/libfftw3.a
energy_LDADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a
energy_LDADD+=@top_builddir@/extlib/Utils/libutilsic.a
