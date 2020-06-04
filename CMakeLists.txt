cmake_minimum_required(VERSION 3.15)
project(sixte)

include_directories($ENV{SIMPUT}/include/)
link_directories($ENV{SIMPUT}/lib/)

set(CMAKE_SIXTE_LOCAL "-g -D_GNU_SOURCE -DHTRS_ARCPIXELS -DSIXT_DATA_PATH='\"$ENV{SIXTE}/share/sixte\"'")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_SIXTE_LOCAL}")

message(" *** Using Simput at $ENV{SIMPUT} ")
message(" *** Using Sixte at  $ENV{SIXTE} ")

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(extlib/progressbar/include)
include_directories(libsixt)
include_directories(tools/athenawfisim)
include_directories(tools/comabackpro)
include_directories(tools/comadet)
include_directories(tools/comaimg)
include_directories(tools/comaimgPM)
include_directories(tools/comaphovign)
include_directories(tools/comarecon)
include_directories(tools/epicmos1_events)
include_directories(tools/epicmos2_events)
include_directories(tools/epicpn_events)
include_directories(tools/ero_calevents)
include_directories(tools/ero_rawevents)
include_directories(tools/erosim)
include_directories(tools/evpat)
include_directories(tools/exposure_map)
include_directories(tools/fudgexp)
include_directories(tools/gendetsim)
include_directories(tools/gennoisespec)
include_directories(tools/gradeddetection)
include_directories(tools/htrssim)
include_directories(tools/ladsim)
include_directories(tools/makelc)
include_directories(tools/makespec)
include_directories(tools/nustarsim)
include_directories(tools/pha2pi)
include_directories(tools/phogen)
include_directories(tools/phoimg)
include_directories(tools/pixdetillum)
include_directories(tools/piximpacts)
include_directories(tools/pulsetemplgen)
include_directories(tools/pulsetemplimport)
include_directories(tools/radec2xy)
include_directories(tools/runmask)
include_directories(tools/runsixt)
include_directories(tools/runtes)
include_directories(tools/streamtotriggers)
include_directories(tools/tes_grades)
include_directories(tools/tesconstpileup)
include_directories(tools/tesgenimpacts)
include_directories(tools/tesreconstruction)
include_directories(tools/tessim)
include_directories(tools/tesstream)
include_directories(tools/xifupipeline)
include_directories(tools/xml2svg)
include_directories(tools/xms_pixtemp)

set(SOURCE_FILES
        libsixt/advdet.c
        libsixt/advdet.h
        libsixt/arcpixels.c
        libsixt/arcpixels.h
        libsixt/attitude.c
        libsixt/attitude.h
        libsixt/attitudefile.c
        libsixt/attitudefile.h
        libsixt/background.c
        libsixt/background.h
        libsixt/badpixmap.c
        libsixt/badpixmap.h
        libsixt/balancing.c
        libsixt/balancing.h
        libsixt/check_fov.c
        libsixt/check_fov.h
        libsixt/clocklist.c
        libsixt/clocklist.h
        libsixt/codedmask.c
        libsixt/codedmask.h
        libsixt/comadetector.c
        libsixt/comadetector.h
        libsixt/comaevent.c
        libsixt/comaevent.h
        libsixt/comaeventfile.c
        libsixt/comaeventfile.h
        libsixt/crosstalk.c
        libsixt/crosstalk.h
        libsixt/det_phi_max.c
        libsixt/det_phi_max.h
        libsixt/detstruct2obj2d.c
        libsixt/detstruct2obj2d.h
        libsixt/event.c
        libsixt/event.h
        libsixt/eventarray.c
        libsixt/eventarray.h
        libsixt/eventfile.c
        libsixt/eventfile.h
        libsixt/eventlist.c
        libsixt/eventlist.h
        libsixt/exponentialchargecloud.h
        libsixt/fft_array.c
        libsixt/fft_array.h
        libsixt/find_position.c
        libsixt/find_position.h
        libsixt/gaussianchargecloud.h
        libsixt/gendet.c
        libsixt/gendet.h
        libsixt/gendetline.c
        libsixt/gendetline.h
        libsixt/genericdetector.c
        libsixt/genericdetector.h
        libsixt/geninst.c
        libsixt/geninst.h
        libsixt/genpixgrid.c
        libsixt/genpixgrid.h
        libsixt/gentel.c
        libsixt/gentel.h
        libsixt/genutils.cpp
        libsixt/genutils.h
        libsixt/grading.c
        libsixt/grading.h
        libsixt/gti.c
        libsixt/gti.h
        libsixt/hexagonalpixels.c
        libsixt/hexagonalpixels.h
        libsixt/htrsdetector.c
        libsixt/htrsdetector.h
        libsixt/htrsevent.h
        libsixt/htrseventfile.c
        libsixt/htrseventfile.h
        libsixt/htrstelstream.c
        libsixt/htrstelstream.h
        libsixt/impact.c
        libsixt/impact.h
        libsixt/impactfile.c
        libsixt/impactfile.h
        libsixt/inoututils.cpp
        libsixt/inoututils.h
        libsixt/integraSIRENA.cpp
        libsixt/integraSIRENA.h
        libsixt/kdtreeelement.c
        libsixt/kdtreeelement.h
        libsixt/lad.c
        libsixt/lad.h
        libsixt/lad_init.c
        libsixt/ladevent.c
        libsixt/ladevent.h
        libsixt/ladeventfile.c
        libsixt/ladeventfile.h
        libsixt/ladimpact.c
        libsixt/ladimpact.h
        libsixt/ladimpactfile.c
        libsixt/ladimpactfile.h
        libsixt/ladsignal.c
        libsixt/ladsignal.h
        libsixt/ladsignalfile.c
        libsixt/ladsignalfile.h
        libsixt/ladsignallist.c
        libsixt/ladsignallist.h
        libsixt/linkedimplist.c
        libsixt/linkedimplist.h
        libsixt/linkedpholist.c
        libsixt/linkedpholist.h
        libsixt/log.cpp
        libsixt/log.h
        libsixt/maskshadow.c
        libsixt/maskshadow.h
        libsixt/masksystem.c
        libsixt/masksystem.h
        libsixt/mt19937ar.c
        libsixt/mt19937ar.h
        libsixt/mxs.c
        libsixt/mxs.h
        libsixt/namelist.c
        libsixt/namelist.h
        libsixt/obj2d.c
        libsixt/obj2d.h
        libsixt/optimalfilters.c
        libsixt/optimalfilters.h
        libsixt/pha2pilib.c
        libsixt/pha2pilib.h
        libsixt/phabkg.c
        libsixt/phabkg.h
        libsixt/phdet.c
        libsixt/phdet.h
        libsixt/phgen.c
        libsixt/phgen.h
        libsixt/phimg.c
        libsixt/phimg.h
        libsixt/photon.c
        libsixt/photon.h
        libsixt/photonfile.c
        libsixt/photonfile.h
        libsixt/phpat.c
        libsixt/phpat.h
        libsixt/phproj.c
        libsixt/phproj.h
        libsixt/pixelimpact.h
        libsixt/pixelimpactfile.c
        libsixt/pixelimpactfile.h
        libsixt/point.h
        libsixt/projectedmask.c
        libsixt/projectedmask.h
        libsixt/psf.c
        libsixt/psf.h
        libsixt/pulseprocess.cpp
        libsixt/pulseprocess.h
        libsixt/radec2xylib.c
        libsixt/radec2xylib.h
        libsixt/reconstruction.c
        libsixt/reconstruction.h
        libsixt/repix.c
        libsixt/repix.h
        libsixt/rndgen.c
        libsixt/rndgen.h
        libsixt/scheduler.cpp
        libsixt/scheduler.h
        libsixt/sixt.c
        libsixt/sixt.h
 #       libsixt/sixt_main.c
        libsixt/sixtesvg.c
        libsixt/sixtesvg.h
        libsixt/skyimage.c
        libsixt/skyimage.h
        libsixt/source.c
        libsixt/source.h
        libsixt/sourcecatalog.c
        libsixt/sourcecatalog.h
        libsixt/sourceimage.c
        libsixt/sourceimage.h
        libsixt/squarepixels.c
        libsixt/squarepixels.h
        libsixt/tasksSIRENA.cpp
        libsixt/tasksSIRENA.h
        libsixt/telemetrypacket.c
        libsixt/telemetrypacket.h
        libsixt/telescope.h
        libsixt/tescrosstalk.c
        libsixt/tescrosstalk.h
        libsixt/tesdatastream.c
        libsixt/tesdatastream.h
        libsixt/teseventlist.c
        libsixt/teseventlist.h
        libsixt/tesinitialization.c
        libsixt/tesinitialization.h
        libsixt/tesnoisespectrum.c
        libsixt/tesnoisespectrum.h
        libsixt/tespixel.h
        libsixt/tesproftemplates.c
        libsixt/tesproftemplates.h
        libsixt/tesreconstruction.h
        libsixt/tesrecord.c
        libsixt/tesrecord.h
        libsixt/tesstreamfile.c
        libsixt/tesstreamfile.h
        libsixt/testimg.c
        libsixt/testimg.h
        libsixt/testrigger.c
        libsixt/testrigger.h
        libsixt/testriggerfile.c
        libsixt/testriggerfile.h
        libsixt/threadsafe_queue.h
        libsixt/vignetting.c
        libsixt/vignetting.h
        libsixt/xmlbuffer.c
        libsixt/xmlbuffer.h
#        test/unit/random_number_gen.c
#        test/unit/test_genpixgrid.c
#        test/unit/test_vignetting.c
#        test/unit/unit_test_all.c
#        config.h
        sixteconfig.h)


#tools/athenawfisim/athenawfisim.c
#tools/athenawfisim/athenawfisim.h
#tools/attgen_dither/attgen_dither.c
#tools/comabackpro/comabackpro.c
#tools/comabackpro/comabackpro.h
#tools/comadet/comadet.c
#tools/comadet/comadet.h
#tools/comaexp/comaexp.c
#tools/comaimg/comaimg.c
#tools/comaimg/comaimg.h
#tools/comaimgPM/comaimgPM.c
#tools/comaimgPM/comaimgPM.h
#tools/comaphovign/comaphovign.c
#tools/comaphovign/comaphovign.h
#tools/comarecon/comarecon.c
#tools/comarecon/comarecon.h
#tools/epicmos1_events/epicmos1_events.c
#tools/epicmos1_events/epicmos1_events.h
#tools/epicmos2_events/epicmos2_events.c
#tools/epicmos2_events/epicmos2_events.h
#tools/epicpn_events/epicpn_events.c
#tools/epicpn_events/epicpn_events.h
#tools/ero_calevents/ero_calevents.c
#tools/ero_calevents/ero_calevents.h
#tools/ero_exposure/ero_exposure.c
#tools/ero_fits2tm/ero_fits2tm.c
#tools/ero_rawevents/ero_rawevents.c
#tools/ero_rawevents/ero_rawevents.h
#tools/ero_vis/ero_vis.c
#tools/erosim/erosim.c
#tools/erosim/erosim.h
#tools/evpat/evpat.c
#tools/evpat/evpat.h
#tools/exposure_map/exposure_map.c
#tools/exposure_map/exposure_map.h
#tools/fudgexp/fudgexp.c
#tools/fudgexp/fudgexp.h
#tools/gendetsim/gendetsim.c
#tools/gendetsim/gendetsim.h
#tools/gennoisespec/gennoisespec.cpp
#tools/gennoisespec/gennoisespec.h
#tools/gradeddetection/gradeddetection.c
#tools/gradeddetection/gradeddetection.h
#tools/htrssim/htrssim.c
#tools/htrssim/htrssim.h
#tools/imgev/imgev.c
#tools/ladsim/ladsim.c
#tools/ladsim/ladsim.h
#tools/makelc/makelc.c
#tools/makelc/makelc.h
#tools/makespec/makespec.c
#tools/makespec/makespec.h
#tools/nustarsim/nustarsim.c
#tools/nustarsim/nustarsim.h
#tools/orbatt/orbatt.c
#tools/pha2pi/pha2pi.c
#tools/pha2pi/pha2pi.h
#tools/phogen/phogen.c
#tools/phogen/phogen.h
#tools/phoimg/phoimg.c
#tools/phoimg/phoimg.h
#tools/pixdetillum/pixdetillum.c
#tools/pixdetillum/pixdetillum.h
#tools/piximpacts/piximpacts.c
#tools/piximpacts/piximpacts.h
#tools/projev/projev.c
#tools/psfgen/psfgen.c
#tools/pulsetemplgen/pulsetemplgen.c
#tools/pulsetemplgen/pulsetemplgen.h
#tools/pulsetemplimport/pulsetemplimport.c
#tools/pulsetemplimport/pulsetemplimport.h
#tools/radec2xy/radec2xy.c
#tools/radec2xy/radec2xy.h
#tools/runmask/runmask.c
#tools/runmask/runmask.h
##        tools/runsixt/runsixt.c
##        tools/runsixt/runsixt.h
#tools/runtes/runtes.c
#tools/runtes/runtes.h
#tools/sixteversion/sixteversion.c
#tools/streamtotriggers/streamtotriggers.c
#tools/streamtotriggers/streamtotriggers.h
#tools/tes_grades/tes_grades.c
#tools/tes_grades/tes_grades.h
#tools/tesconstpileup/tesconstpileup.c
#tools/tesconstpileup/tesconstpileup.h
#tools/tesgenimpacts/tesgenimpacts.c
#tools/tesgenimpacts/tesgenimpacts.h
#tools/tesreconstruction/tesreconstruction.c
#tools/tesreconstruction/tesreconstruction.h
#tools/tessim/tes_models.c
#tools/tessim/tes_simulation.c
#tools/tessim/tessim.c
#tools/tessim/tessim.h
#tools/tessim/tessim_bbfb.c
#tools/tessim/tessim_datastream.c
#tools/tessim/tessim_impactlist.c
#tools/tessim/tessim_solvers.c
#tools/tessim/tessim_tesrecord.c
#tools/tessim/tessim_trigger.c
#tools/tesstream/tesstream.c
#tools/tesstream/tesstream.h
#tools/xifupipeline/xifupipeline.c
#tools/xifupipeline/xifupipeline.h
#tools/xml2svg/xml2svg.c
#tools/xml2svg/xml2svg.h
#tools/xms_pixtemp/xms_pixtemp.c
#tools/xms_pixtemp/xms_pixtemp.h


set(EXEC_FILES runsixt )

set(TOOL_DIR tools)


set(SIMPUT_LIBS cfitsio simput ape atFunctions fftw3 hdinit hdio hdsp hdutils labnh posstring wcs)
set(EXT_LIBS expat readline ncurses m gsl gslcblas)

# special libraries
find_package (Threads)
set(EXT_LIBS ${EXT_LIBS} ${CMAKE_THREAD_LIBS_INIT})

foreach (execfile ${EXEC_FILES})
    message(${execfile})
    add_executable(${execfile} ${TOOL_DIR}/${execfile}/${execfile}.c ${TOOL_DIR}/${execfile}/${execfile}.h ${SOURCE_FILES})
    target_link_libraries(${execfile} ${EXT_LIBS} ${SIMPUT_LIBS} )
endforeach (execfile ${EXEC_FILES})
