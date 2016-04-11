#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/srg

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml'"  > /dev/stderr    
    echo "  --> skip testing 'erosim', 'ero_calevents',and 'ero_rawevents' "  > /dev/stderr    
    exit 
endif

erosim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/erosim_ \
  XMLFile1=${xml}/erosita_1.xml \
  XMLFile2=${xml}/erosita_2.xml \
  XMLFile3=${xml}/erosita_3.xml \
  XMLFile4=${xml}/erosita_4.xml \
  XMLFile5=${xml}/erosita_5.xml \
  XMLFile6=${xml}/erosita_6.xml \
  XMLFile7=${xml}/erosita_7.xml \
  Simput=$indir/simput.fits \
  Exposure=1 \
  clobber=yes 

set num_wc = `ls output/erosim*_evt.fits* | wc -l`
if ($num_wc != 7) then
    echo " *** error in erosim *** not the expected number of files created"  > /dev/stderr
    set status = 1
    exit
endif


erosim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/keep_erosim_ \
  RawData=raw.fits \
  XMLFile1=${xml}/erosita_1.xml \
  XMLFile2=${xml}/erosita_2.xml \
  XMLFile3=${xml}/erosita_3.xml \
  XMLFile4=${xml}/erosita_4.xml \
  XMLFile5=${xml}/erosita_5.xml \
  XMLFile6=${xml}/erosita_6.xml \
  XMLFile7=${xml}/erosita_7.xml \
  Simput=$indir/simput.fits \
  Exposure=1 \
  clobber=yes 

set num_wc = `ls output/keep_erosim*.fits | wc -l`
if ($num_wc != 14) then
    echo " *** error in erosim *** not the expected number of files created"  > /dev/stderr
    set status = 1
    exit
endif


ero_calevents \
  EvtFile=$outdir/erosim_ccd1_evt.fits \
  eroEvtFile=$outdir/erosim_ccd1_eroevt.fits \
  CCDNr=1 \
  RA=0.0 Dec=0.0 \
  clobber=yes

set num_wc = `ls output/erosim_ccd1_eroevt.fits | wc -l`
if ($num_wc != 1) then
    echo " *** error in erosim *** not the expected number of files created"  > /dev/stderr
    set status = 1
    exit
endif


ero_rawevents \
  RawData=$outdir/keep_erosim_ccd1_raw.fits \
  eroEvtFile=$outdir/erosim_ccd1_raweroevt.fits \
  CCDNr=1 \
  clobber=yes

set num_wc = `ls output/erosim_ccd1_raweroevt.fits | wc -l`
if ($num_wc != 1) then
    echo " *** error in erosim *** not the expected number of files created"  > /dev/stderr
    set status = 1
    exit
endif
