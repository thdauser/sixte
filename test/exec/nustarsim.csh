#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/nustar/

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'nustarsim' "      
    exit 
endif

#gdb --args  
nustarsim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/nustarsim_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/nustar.xml \
  Exposure=1 \
  clobber=yes 


nustarsim  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/keep_nustarsim_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/nustar.xml \
  Exposure=1 \
  clobber=yes 
