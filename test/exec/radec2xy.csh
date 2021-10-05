#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena-wfi/wfi_wo_filter_baseline

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'makespec' "
    exit
endif

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/radec2xy_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/ld_wfi_ff_large.xml \
  Exposure=100 \
  clobber=yes

radec2xy \
  EvtFile=$outdir/radec2xy_evt.fits \
  RefRa=0. \
  RefDec=0. \
  Projection="AIT"

radec2xy \
  EvtFile=$outdir/radec2xy_evt.fits \
  RefRa=10. \
  RefDec=10. \
  Projection="SIN"
