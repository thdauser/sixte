AM_CFLAGS=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

if RCL
AM_CFLAGS+=-I$(RCL_INC) -DUSE_RCL
AM_LDFLAGS=-L$(RCL_LIB) -lrcl
endif

########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir=$(pkgdatadir)/templates
dist_fitstemplates_DATA =coma.eventlist.tpl
dist_fitstemplates_DATA+=photonfile.tpl 
dist_fitstemplates_DATA+=impactfile.tpl
dist_fitstemplates_DATA+=eventfile.tpl
dist_fitstemplates_DATA+=ladimpactfile.tpl
dist_fitstemplates_DATA+=ladsignalfile.tpl
dist_fitstemplates_DATA+=ladeventfile.tpl
dist_fitstemplates_DATA+=makespec.tpl
dist_fitstemplates_DATA+=makelc.tpl
dist_fitstemplates_DATA+=htrs.eventlist.tpl
dist_fitstemplates_DATA+=pixelimpactfile.tpl

# Directory where to install the XML files.
xmldir=$(pkgdatadir)/xml
dist_xml_DATA =wfi.iaat.xml
dist_xml_DATA+=htrs.xml
dist_xml_DATA+=comadet.xml
dist_xml_DATA+=aciss3.xml


########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libsixt.la
libsixt_la_LDFLAGS = -version-info 2:0:0

# Sources:
libsixt_la_SOURCES=genericdetector.c squarepixels.c eventlist.c		\
		  impactfile.c ladimpactfile.c htrsdetector.c		\
		  geninst.c gendet.c gentel.c genpixgrid.c		\
		  gendetline.c ladsignalfile.c ladeventfile.c		\
		  phabkg.c eventfile.c clocklist.c badpixmap.c		\
		  htrseventfile.c hexagonalpixels.c arcpixels.c		\
		  telemetrypacket.c htrstelstream.c comadetector.c	\
		  comaeventfile.c psf.c vignetting.c codedmask.c	\
		  attitude.c attitudefile.c sixt.c photon.c		\
		  check_fov.c photonfile.c kdtreeelement.c		\
		  sourcecatalog.c source.c linkedpholist.c		\
		  ladsignallist.c erodetbkgrndgen.c phgen.c phimg.c	\
		  phdet.c phproj.c phpat.c event.c ladsignal.c		\
		  ladevent.c ladimpact.c lad.c lad_init.c xmlbuffer.c	\
		  gti.c sourceimage.c reconstruction.c eventarray.c	\
		  fft_array.c balancing.c find_position.c det_phi_max.c \
	          advdet.c pixelimpactfile.c tesdatastream.c 		\
		  tesnoisespectrum.c tesproftemplates.c maskshadow.c 	\
		  tesstreamfile.c testriggerfile.c testimg.c 		\
		  projectedmask.c repix.c tesinitialization.c	        \
		  comaevent.c skyimage.c detstruct2obj2d.c obj2d.c 	\
		  sixtesvg.c tesrecord.c teseventlist.c
		  
############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=genericdetector.h squarepixels.h hexagonalpixels.h	\
		gaussianchargecloud.h exponentialchargecloud.h		\
		arcpixels.h eventlist.h impactfile.h ladimpactfile.h	\
		impact.h geninst.h gendet.h gentel.h genpixgrid.h	\
		gendetline.h clocklist.h ladsignal.h ladevent.h		\
		ladimpact.h event.h ladsignalfile.h phabkg.h		\
		ladeventfile.h eventfile.h badpixmap.h htrsdetector.h	\
		htrseventfile.h htrsevent.h telemetrypacket.h		\
		htrstelstream.h comadetector.h comaeventfile.h		\
		comaevent.h psf.h vignetting.h codedmask.h attitude.h	\
		attitudefile.h telescope.h sixt.h point.h photon.h	\
		check_fov.h photonfile.h kdtreeelement.h		\
		sourcecatalog.h source.h linkedpholist.h		\
		ladsignallist.h erodetbkgrndgen.h phgen.h phimg.h	\
		phdet.h phproj.h phpat.h lad.h xmlbuffer.h gti.h	\
		sourceimage.h reconstruction.h eventarray.h		\
		fft_array.h balancing.h find_position.h			\
		det_phi_max.h advdet.h 					\
		pixelimpactfile.h pixelimpact.h tesdatastream.h 	\
		tesnoisespectrum.h tesproftemplates.h maskshadow.h      \
		tesstreamfile.h testriggerfile.h testimg.h  		\
		projectedmask.h repix.h tesinitialization.h skyimage.h	\
		detstruct2obj2d.h obj2d.h sixtesvg.h tesrecord.h	\
		teseventlist.h
