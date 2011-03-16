AM_CFLAGS  = -I${HEADAS}/include
#AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = coma.eventlist.tpl
dist_fitstemplates_DATA += erosita.eventlist.tpl 
dist_fitstemplates_DATA += eventlistfile.tpl
dist_fitstemplates_DATA += genpatternfile.tpl
dist_fitstemplates_DATA += htrs.eventlist.tpl
dist_fitstemplates_DATA += impactlist.tpl
dist_fitstemplates_DATA += photonlist.tpl 
dist_fitstemplates_DATA	+= wfi.full1024.eventlist.tpl 
dist_fitstemplates_DATA += wfi.window16.eventlist.tpl
dist_fitstemplates_DATA += wfi.labor64.eventlist.tpl
dist_fitstemplates_DATA += xms.eventlist.tpl

# Directory where to install the XML files.
xmldir = $(pkgdatadir)/xml
dist_xml_DATA  = eroframestore.xml 
dist_xml_DATA += wfi.full.xml
dist_xml_DATA += wfi.iaat.xml
dist_xml_DATA += htrs.xml
dist_xml_DATA += epicpn.fullframe.xml
dist_xml_DATA += epicpn.timing.xml
dist_xml_DATA += epicpn.burst.xml
dist_xml_DATA += comadet.xml
dist_xml_DATA += gravitas.xml
dist_xml_DATA += aciss3.xml

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
		     eventlistfile.c genpatternfile.c clocklist.c	\
		     badpixmap.c geneventgrading.c htrseventfile.c	\
		     hexagonalpixels.c arcpixels.c telemetrypacket.c	\
		     htrstelstream.c comadetector.c comaeventfile.c	\
		     epicpn.c psf.c vignetting.c codedmask.c		\
		     attitudecatalog.c attitudefile.c sixt.c vector.c	\
		     pointsources.c pointsourcefile.c			\
		     pointsourcelist.c pointsourcecatalog.c		\
		     extendedsources.c sourceimage.c photon.c		\
		     spectrum.c fits_pha.c check_fov.c			\
		     photonlistfile.c linlightcurve.c kdtree.c rmf.c	\
		     arf.c kdtreeelement.c xraysourcecatalog.c		\
		     xraysourcespectrum.c xraysource.c linkedpholist.c	\
		     erodetbkgrndgen.c phgen.c phimg.c phdet.c


############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS = genericdetector.h squarepixels.h hexagonalpixels.h	\
		  gaussianchargecloud.h exponentialchargecloud.h	\
		  arcpixels.h eventfile.h impactlistfile.h impact.h	\
		  uniformdetectorbackground.h framestoredetector.h	\
		  erositaeventfile.h erositaevent.h wfidetector.h	\
		  wfieventfile.h wfievent.h xmsdetector.h gendet.h	\
		  genpixgrid.h gensplit.h gendetline.h clocklist.h	\
		  event.h genpattern.h eventlistfile.h genpatternfile.h	\
		  genpileupflag.h badpixmap.h geneventgrading.h		\
		  xmseventfile.h xmsevent.h htrsdetector.h		\
		  htrseventfile.h htrsevent.h telemetrypacket.h		\
		  htrstelstream.h comadetector.h comaeventfile.h	\
		  comaevent.h epicpn.h psf.h vignetting.h codedmask.h	\
		  attitudecatalog.h attitudefile.h telescope.h sixt.h	\
		  point.h vector.h pointsources.h pointsourcefile.h	\
		  pointsourcelist.h pointsourcecatalog.h		\
		  extendedsources.h sourceimage.h photon.h spectrum.h	\
		  fits_pha.h check_fov.h photonlistfile.h		\
		  linlightcurve.h kdtree.h rmf.h arf.h kdtreeelement.h	\
		  xraysourcecatalog.h xraysourcespectrum.h		\
		  xraysource.h linkedpholist.h erodetbkgrndgen.h	\
		  phgen.h phimg.h phdet.h

