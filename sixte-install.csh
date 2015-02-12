#!/bin/csh

# install for SIXTE - csh version
#
# This script assumes that the SIXTE environment variable is set
# and uses nifty tricks from HEADAS' setup scripts
#
# Author: Joern Wilms, joern.wilms@sternwarte.uni-erlangen.de
#

if (${?SIXTE} == 0) then
  echo "sixte-install.csh: ERROR -- set SIXTE before sourcing sixte-install.csh"
  exit 1
endif

if (! -d ${SIXTE}) then
    echo "Directory ${SIXTE} does not exist"
    exit 2
endif

#
# set the SIMPUT environment variable
#

if (${?SIMPUT} == 0) then
  setenv SIMPUT ${SIXTE}
    if (! -e ${SIMPUT}/bin/simput-install.csh) then
	echo "sixte-install.csh: ERROR -- set SIMPUT environment variable before sourcing sixte-install.csh"
	exit 1
    endif
endif

#
# run the setup script for SIMPUT
#
source ${SIMPUT}/bin/simput-install.csh

#
# set paths
#
set SIXTE_BIN = ${SIXTE}/bin
setenv PATH ${SIXTE_BIN}:${PATH}

#
# setup parameter files
#
if (${?PFILES} == 0) then
    mkdir -p ${HOME}/pfiles
    setenv PFILES "${HOME}/pfiles;${SIXTE}/share/sixte/pfiles"
else
    setenv PFILES "${PFILES}:${SIXTE}/share/sixte/pfiles"
endif

#
# set LD_LIBRARY_PATH
#
set SIMPUT_LIB = ${SIMPUT}/lib

if (${?LD_LIBRARY_PATH} == 0) then
    setenv LD_LIBRARY_PATH ${SIMPUT_LIB}
else
    setenv LD_LIBRARY_PATH `echo ":${LD_LIBRARY_PATH}:" | sed "s%:${SIMPUT_LIB}:%:%g" | sed 's%::*$%%'`
    setenv LD_LIBRARY_PATH ${SIMPUT_LIB}${LD_LIBRARY_PATH}
endif

set build_os = `uname`
if (${build_os} == "Darwin") then
    if (${?DYLD_LIBRARY_PATH} == 0) then
	setenv DYLD_LIBRARY_PATH ${SIMPUT_LIB}
    else
	setenv DYLD_LIBRARY_PATH `echo ":${DYLD_LIBRARY_PATH}:" | sed "s%:${SIMPUT_LIB}:%:%g" | sed 's%::*$%%'`
	setenv DYLD_LIBRARY_PATH ${SIMPUT_LIB}${DYLD_LIBRARY_PATH}
    endif
endif
