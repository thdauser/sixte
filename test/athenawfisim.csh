#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena-wfi/wfi_w_filter

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'athenawfisim' "      
    exit 
endif

athenawfisim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/athenawfisim_ \
  EvtFile=evt.fits \
  Simput=$indir/simput.fits \
  XMLFile0=${xml}/ld_wfi_ff_chip0.xml \
  XMLFile1=${xml}/ld_wfi_ff_chip1.xml \
  XMLFile2=${xml}/ld_wfi_ff_chip2.xml \
  XMLFile3=${xml}/ld_wfi_ff_chip3.xml \
  Exposure=1 \
  clobber=yes 


athenawfisim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/keep_athenawfisim_ \
  EvtFile=evt.fits \
  RawData=raw.fits \
  Simput=$indir/simput.fits \
  XMLFile0=${xml}/ld_wfi_ff_chip0.xml \
  XMLFile1=${xml}/ld_wfi_ff_chip1.xml \
  XMLFile2=${xml}/ld_wfi_ff_chip2.xml \
  XMLFile3=${xml}/ld_wfi_ff_chip3.xml \
  Exposure=1 \
  clobber=yes 
