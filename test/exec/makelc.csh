#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena/1469mm_wfi_wo_filter

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'makelc' "
    exit
endif

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/makelc_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/depfet_b_1l_ff_large.xml \
  Exposure=100 \
  clobber=yes

makelc \
  EvtFile=$outdir/makelc_evt.fits \
  LightCurve=$outdir/makelc_lc.fits \
  length=100.0 \
  dt=1.0 \
  clobber=yes
