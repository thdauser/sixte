#!/bin/bash -v

#################################################################
# Correcting the dependencies in SIXTE libraries
# Author: E. Cucchetti (IRAP) & Thomas Dauser (ECAP), 2017-11-13
# If problems remains, contact sixte-support@lists.fau.de
#################################################################

install_prefix_simput=$1
install_prefix_sixte=$2

for tool in bin/sixte_arfgen bin/athenawfisim bin/comabackpro bin/comadet bin/comaexp bin/comaimg bin/comaimgPM \
    bin/comaphovign bin/comarecon bin/epicmos1_events bin/epicmos2_events bin/ero_calevents \
    bin/ero_exposure bin/ero_fits2tm bin/ero_rawevents bin/ero_vis bin/erosim bin/evpat \
    bin/exposure_map bin/fudgexp bin/gendetsim bin/gennoisespec bin/gradeddetection \
    bin/htrssim bin/imgev bin/ladsim bin/makelc bin/makespec bin/nustarsim bin/orbatt \
    bin/phogen bin/phoimg bin/pixdetillum bin/piximpacts bin/projev bin/psfgen bin/pulsetemplgen \
    bin/pulsetemplimport bin/runsixt bin/runtes bin/streamtotriggers bin/tes_grades bin/tesconstpileup \
    bin/tesgenimpacts bin/tesreconstruction bin/tessim bin/tesstream bin/xifupipeline \
    bin/xml2svg bin/xms_pixtemp lib/libsixt.2.dylib lib/libprogressbar.0.dylib
do
    install_name_tool -change @rpath/libcfitsio.9.dylib $install_prefix_simput/lib/libcfitsio.9.dylib $install_prefix_sixte/$tool
    install_name_tool -change libwcs.7.dylib $install_prefix_simput/lib/libwcs.7.dylib $install_prefix_sixte/$tool
done
