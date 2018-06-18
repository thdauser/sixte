#!/bin/csh

source setup/setup.csh

set xml = $xmldir/athena-xifu
set std = $xml/xifu_baseline.xml
set adv = $xml/xifu_detector_lpa2_pixoffset_mux40.xml

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml'"      
    echo "  --> skip testing 'xifupipeline' "      
    exit 
endif

xifupipeline  \
    Prefix=$outdir/xifupipeline_ \
    XMLFile=$std \
    AdvXml=$adv \
    RA=0.0 Dec=0.0 \
    Background=no \
    Simput=$indir/simput.fits \
    Exposure=10. \
    UseRMF=yes \
    doCrosstalk=none \
    chatter=3 \
    clobber=yes


xifupipeline  \
    Prefix=$outdir/keep_xifupipeline_ \
    PixImpactList=piximp_raw.fits \
    XMLFile=$std \
    AdvXml=$adv \
    RA=0.0 Dec=0.0 \
    Background=no \
    Simput=$indir/simput.fits \
    Exposure=10. \
    UseRMF=yes \
    doCrosstalk=no \
    chatter=3 \
    clobber=yes

    
