#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/athena-wfi/wfi_wo_filter_baseline

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' : skip testing 'imgev' "
    exit
endif

runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/imgev_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/ld_wfi_ff_large.xml \
  Exposure=100 \
  clobber=yes

imgev \
  EvtFile=$outdir/imgev_evt.fits \
  Image=$outdir/imgev_img.fits \
  CoordinateSystem=0 \
  Projection=AIT \
  NAXIS1=36 NAXIS2=18 \
  CUNIT1=deg CUNIT2=deg \
  CRPIX1=18.5 CRPIX2=9.5 \
  CRVAL1=0.0 CRVAL2=0.0 \
  CDELT1=-0.5 CDELT2=0.5 \
  clobber=yes
