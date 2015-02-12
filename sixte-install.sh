#!/bin/sh

# install for SIXTE - sh version
#
# This script assumes that the SIXTE environment variable is set
# and uses nifty tricks from HEADAS' setup scripts
#
# Author: Joern Wilms, joern.wilms@sternwarte.uni-erlangen.de
#

if [ "X${SIXTE}" == X ];  then
  echo "sixte-install.csh: ERROR -- set SIXTE before sourcing sixte-install.csh"
  exit 1
fi

if ! [ -d ${SIXTE} ]; then
    echo "Directory ${SIXTE} does not exist"
    exit 2
fi

if [ "X${SIMPUT}" == X ];  then
    export SIMPUT=${SIXTE}
    if ! [ -e ${SIMPUT}/bin/simput-install.csh ]; then
	echo "sixte-install.csh: ERROR -- set SIMPUT environment variable before sourcing sixte-install.csh"
	exit 1
    fi
fi
	
#
# run the setup script for SIMPUT
#
. ${SIMPUT}/bin/simput-install.sh

#
# set paths
#
sixte_bin=${SIXTE}/bin
PATH=${sixte_bin}:${PATH}

#
# setup parameter files
#
if [ "X${PFILES}" == X ]; then
    mkdir -p ${HOME}/pfiles
    PFILES="${HOME}/pfiles;${SIXTE}/share/sixte/pfiles"
else
    PFILES="${PFILES}:${SIXTE}/share/sixte/pfiles"
fi


export DYLD_LIBRARY_PATH LD_LIBRARY_PATH PATH PFILES

SIXTE_LIB=${SIXTE}/lib
if [ "x$LD_LIBRARY_PATH" = x ]; then
  LD_LIBRARY_PATH="$SIXTE_LIB"
else
  LD_LIBRARY_PATH=`echo ":$LD_LIBRARY_PATH:" | sed "s%:$SIXTE_LIB:%:%g" | sed "s%::*$%%"`
  LD_LIBRARY_PATH="$SIXTE_LIB$LD_LIBRARY_PATH"
fi

build_os=`uname`
if [ "$build_os" == "Darwin" ]; then
    if [ "x$DYLD_LIBRARY_PATH" = x ]; then
	DYLD_LIBRARY_PATH="$SIXTE_LIB"
    else
	DYLD_LIBRARY_PATH=`echo ":$DYLD_LIBRARY_PATH:" | sed "s%:$SIXTE_LIB:%:%g" | sed "s%::*$%%"`
	DYLD_LIBRARY_PATH="$SIXTE_LIB$DYLD_LIBRARY_PATH"
    fi
fi
