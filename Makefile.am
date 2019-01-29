dist_bin_SCRIPTS =sixte-install.sh
dist_bin_SCRIPTS+=sixte-install.csh

ACLOCAL_AMFLAGS = -I build-aux/m4

EXTRA_DIST=INSTALL.txt mac_solve_sixte_dependencies.sh # include INSTALL.txt in the distributed package

# The sub-directories are built before the current directory.
# In order to change this, include "." in the list of SUBDIRS.
SUBDIRS=extlib libsixt tools test

# need also to tell distcheck where simput is
AM_DISTCHECK_CONFIGURE_FLAGS= --with-simput=${SIMPUT}

# Define the headers that will be installed in $(includedir):
include_HEADERS=sixteconfig.h

.PHONY: test
test:
	cd test/exec/ && make test
	cd test/e2e/ && ./run_e2e_test.sh

crit_mac_version=9
install-exec-hook:
	@if [ 0$(OSX_VERSION_MINOR) -gt $(crit_mac_version)  ]; then\
		echo "\n  ***  Mac Version $(OSX_VERSION_MAJOR).$(OSX_VERSION_MINOR) detected *** ";\
		echo "\n   ... trying to solve known problems with this installation with this script:\n\n ";\
		./mac_solve_sixte_dependencies.sh @simput_dir@ @prefix@;\
	fi

EXTRA_DIST+=.version
BUILT_SOURCES = $(top_srcdir)/.version
$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@
dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version                                                                       
