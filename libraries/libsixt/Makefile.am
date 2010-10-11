AM_CFLAGS  = -I${HEADAS}/include
#AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = erosita.eventlist.tpl 
dist_fitstemplates_DATA	+= wfi.full1024.eventlist.tpl 
dist_fitstemplates_DATA += wfi.window16.eventlist.tpl
dist_fitstemplates_DATA += wfi.labor64.eventlist.tpl
dist_fitstemplates_DATA += xms.eventlist.tpl
dist_fitstemplates_DATA += htrs.eventlist.tpl
dist_fitstemplates_DATA += impactlist.tpl
dist_fitstemplates_DATA += coma.eventlist.tpl
dist_fitstemplates_DATA += geneventfile.tpl
dist_fitstemplates_DATA += photonlist.tpl 

# Directory where to install the XML files.
xmldir = $(pkgdatadir)/xml
dist_xml_DATA  = eroframestore.xml 
dist_xml_DATA += wfi.full.xml
dist_xml_DATA += wfi.iaat.xml
dist_xml_DATA += epicpn.burst.xml
dist_xml_DATA += epicpn.timing.xml
dist_xml_DATA += comadet.xml

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libsixt.la

# Sources:
libsixt_la_SOURCES = genericdetector.c squarepixels.c eventfile.c	\
		     impactlistfile.c uniformdetectorbackground.c	\
		     framestoredetector.c erositaeventfile.c		\
		     wfidetector.c wfieventfile.c xmsdetector.c		\
		     xmseventfile.c htrsdetector.c gendet.c		\
		     genpixgrid.c gensplit.c gendetline.c		\
		     geneventfile.c clocklist.c htrseventfile.c		\
		     hexagonalpixels.c arcpixels.c telemetrypacket.c	\
		     htrstelstream.c comadetector.c comaeventfile.c	\
		     epicpn.c psf.c vignetting.c codedmask.c		\
		     attitudecatalog.c attitudefile.c sixt_random.c	\
		     sixt_string.c vector.c pointsources.c		\
		     pointsourcefile.c pointsourcelist.c		\
		     pointsourcecatalog.c extendedsources.c		\
		     sourceimage.c photon.c spectrum.c fits_pha.c	\
		     check_fov.c photonlistfile.c linlightcurve.c	\
		     kdtree.c rmf.c arf.c


############ HEADERS #################

# Define the headers that will be installed in $(includedir) :
include_HEADERS = genericdetector.h squarepixels.h hexagonalpixels.h	\
		  gaussianchargecloud.h exponentialchargecloud.h	\
		  arcpixels.h eventfile.h impactlistfile.h impact.h	\
		  uniformdetectorbackground.h framestoredetector.h	\
		  erositaeventfile.h erositaevent.h wfidetector.h	\
		  wfieventfile.h wfievent.h xmsdetector.h gendet.h	\
		  genpixgrid.h gensplit.h gendetline.h clocklist.h	\
		  genevent.h geneventfile.h genpileupflag.h		\
		  xmseventfile.h xmsevent.h htrsdetector.h		\
		  htrseventfile.h htrsevent.h telemetrypacket.h		\
		  htrstelstream.h comadetector.h comaeventfile.h	\
		  comaevent.h epicpn.h psf.h vignetting.h codedmask.h	\
		  attitudecatalog.h attitudefile.h telescope.h sixt.h	\
		  sixt_random.h sixt_string.h point.h vector.h		\
		  astrosources.h pointsources.h pointsourcefile.h	\
		  pointsourcelist.c pointsourcecatalog.h		\
		  extendedsources.h sourceimage.h photon.h spectrum.h	\
		  fits_pha.h check_fov.h photonlistfile.h		\
		  linlightcurve.h kdtree.h rmf.h arf.h

