#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena/1469mm_wfi_wo_filter

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'makespec' "  > /dev/stderr    
    exit 
endif

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/makespec_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/depfet_b_1l_ff_large.xml \
  Mission=Athena \
  Instrument=WFI \
  Mode=large \
  Exposure=100 \
  clobber=yes 

makespec \
  EvtFile=$outdir/makespec_evt.fits \
  Spectrum=$outdir/makespec_pha.fits \
  RSPPath=$xml \
  clobber=yes

