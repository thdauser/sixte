AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/Utils

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
#pfilesdir=$(pkgdatadir)/pfiles
#dist_pfiles_DATA=trigger.par

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
bin_PROGRAMS=trigger

trigger_SOURCES=trigger.cpp trigger.h
trigger_LDADD =@top_builddir@/libsixt/libsixt.la
trigger_LDADD+=@top_builddir@/extlib/Utils/libutilsic.la
