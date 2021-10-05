#!/bin/csh

source setup/setup.csh

set xml = $xmldir/athena-xifu
set adv = $xml/xifu_detector_lpa_75um_AR0.5_pixoffset_mux40_pitch275um.xml

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml'"
    echo "  --> skip testing 'gradeddetection' "
    exit
endif

gradeddetection  \
    PixImpList=$indir/fakeimpact.fits \
    EvtFile=$outdir/evt.fits \
    AdvXml=$adv \
    chatter=3 \
    clobber=yes
