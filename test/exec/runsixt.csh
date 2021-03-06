#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena/1469mm_wfi_wo_filter

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'runsixt' "      
    exit 
endif


#gdb --args  
runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/runsixt_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/depfet_b_1l_ff_large.xml \
  Mission=Athena \
  Instrument=WFI \
  Mode=large \
  Exposure=1 \
  clobber=yes 

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/keep_runsixt_ \
  RawData=raw.fits \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/depfet_b_1l_ff_large.xml \
  Mission=Athena \
  Instrument=WFI \
  Mode=large \
  Exposure=1 \
  clobber=yes 
