#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena-wfi/wfi_wo_filter_baseline

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'runsixt' "
    exit
endif


#gdb --args
runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/runsixt_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/ld_wfi_ff_large.xml \
  Exposure=1 \
  clobber=yes

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/keep_runsixt_ \
  RawData=raw.fits \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/ld_wfi_ff_large.xml \
  Exposure=1 \
  clobber=yes
