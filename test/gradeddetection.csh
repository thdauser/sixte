#!/bin/csh

source setup/setup.csh

set xml = $xmldir/athena/1469mm_xifu
set adv = $xml/xifu_detector_hex_baseline.xml

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml'"  > /dev/stderr    
    echo "  --> skip testing 'gradeddetection' "  > /dev/stderr    
    exit 
endif

gradeddetection  \
    PixImpList=$indir/fakeimpact.fits \
    EvtFile=$outdir/evt.fits \
    AdvXml=$adv \
    chatter=3 \
    clobber=yes
