dist_bin_SCRIPTS =sixte-install.sh
dist_bin_SCRIPTS+=sixte-install.csh

ACLOCAL_AMFLAGS = -I ./build-aux/m4/

EXTRA_DIST=INSTALL.txt # include INSTALL.txt in the distributed package

# The sub-directories are built before the current directory.
# In order to change this, include "." in the list of SUBDIRS.
SUBDIRS=extlib libsixt tools

# need also to tell distcheck where simput is
AM_DISTCHECK_CONFIGURE_FLAGS= --with-libsimput-prefix=${SIMPUT}

# Define the headers that will be installed in $(includedir):
include_HEADERS=sixteconfig.h