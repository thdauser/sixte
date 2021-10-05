#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena-wfi/wfi_wo_filter_baseline

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'makespec' "
    exit
endif

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/makespec_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/ld_wfi_ff_large.xml \
  Exposure=100 \
  clobber=yes

makespec \
  EvtFile=$outdir/makespec_evt.fits \
  Spectrum=$outdir/makespec_pha.fits \
  RSPPath=$xml \
  clobber=yes
