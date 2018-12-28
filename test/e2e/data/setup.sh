#! /bin/bash

. ${SIXTE}/bin/sixte-install.sh

instdir=${SIXTE}/share/sixte/instruments

fname_pholist="pho.fits"
fname_implist="imp.fits"
fname_rawlist="raw.fits"
fname_evtlist="evt.fits"

def_inst="dummy_inst"
def_mode="dummy_mode"
def_filt="dummy_filt"
def_miss="dummy_miss"

def_xml="default_inst.xml"

defRA=0.0
defDec=0.0

defExpos=100

defSimput="dummy.simput" 

defSeed=0

prefix_refdata="ref_"
prefix_dummy="dummy_"
