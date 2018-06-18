#! /bin/csh

source setup/setup.csh
set xml = ${xmldir}/xmm

if (! (-d $xml)) then
    echo " *** warning *** did not find required instrument-directory '$xml' "      
    echo "  --> skip testing 'epicpn_events', 'epicmos1_events', and 'epicmos2_events' "      
    exit 
endif


runsixt  \
  RA=0.0 Dec=0.0 \
  Prefix=$outdir/fakexmm_ \
  Simput=$indir/simput.fits \
  XMLFile=${xml}/epicpn/fullframe_thickfilter.xml \
  Mission=XMM \
  Instrument=EPICpn \
  Mode=ff \
  Exposure=1 \
  MJDREF=50814.0 \
  clobber=yes 

epicpn_events \
  EvtFile=$outdir/fakexmm_evt.fits \
  EPICpnEventList=$outdir/fakexmm_epicpn_evt.fits \
  clobber=yes 

epicmos1_events \
  EvtFile=$outdir/fakexmm_evt.fits \
  EPICmos1EventList=$outdir/fakexmm_epicmos1_evt.fits \
  clobber=yes 

epicmos2_events \
  EvtFile=$outdir/fakexmm_evt.fits \
  EPICmos2EventList=$outdir/fakexmm_epicmos2_evt.fits \
  clobber=yes 

