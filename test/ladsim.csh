#! /bin/csh -fe

source setup/setup.csh
set xml = ${xmldir}/loft/

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'ladsim' "  > /dev/stderr    
    exit 
endif

ladsim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/ladsim_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/lad.xml \
  Exposure=1 \
  clobber=yes 

set num_wc = `ls output/ladsim* | wc -l`
if ($num_wc != 1) then
    echo " *** error in ladsim *** not the expected number of files created"  > /dev/stderr
    set status = 1
    exit
endif

ladsim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/keep_ladsim_ \
  Simput=$indir/simput.fits \
  RawData=raw.fits \
  XMLFile=${xml}/lad.xml \
  Exposure=1 \
  clobber=yes 

set num_wc = `ls output/keep_ladsim* | wc -l`
if ($num_wc != 2) then
    echo " *** error in ladsim *** not the expected number of files created"  > /dev/stderr
    set status = 1
    exit
endif
