#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/srg

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'projev' and 'evpat' "
    exit
endif

runsixt  \
  RA=0.0 Dec=0.0 \
  ImpactList=impact.fits \
  RawData=raw.fits \
  EvtFile=dummy_evt.fits \
  Prefix=$outdir/projev_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/erosita_dummy.xml \
  Exposure=1 \
  clobber=yes

projev  \
  RawData=$outdir/projev_raw.fits \
  RA=0.0 Dec=0.0 \
  XMLFile=${xml}/erosita_dummy.xml \
  Exposure=1 \
  clobber=yes

evpat  \
  RawData=$outdir/projev_raw.fits \
  EvtFile=$outdir/projev_evt.fits \
  XMLFile=${xml}/erosita_dummy.xml \
  clobber=yes
