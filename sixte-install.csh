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
