AM_CFLAGS  = -I${HEADAS}/include
#AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA  = coma.eventlist.tpl
dist_fitstemplates_DATA += eventlist.tpl
dist_fitstemplates_DATA += eroeventlist.tpl
dist_fitstemplates_DATA += patternlist.tpl
dist_fitstemplates_DATA += htrs.eventlist.tpl
dist_fitstemplates_DATA += impactlist.tpl
dist_fitstemplates_DATA += photonlist.tpl 

# Directory where to install the XML files.
xmldir = $(pkgdatadir)/xml
dist_xml_DATA  = wfi.iaat.xml
dist_xml_DATA += htrs.xml
dist_xml_DATA += epicpn.fullframe.xml
dist_xml_DATA += epicpn.timing.xml
dist_xml_DATA += epicpn.burst.xml
dist_xml_DATA += comadet.xml
dist_xml_DATA += aciss3.xml

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libsixt.la

# Sources:
libsixt_la_SOURCES = genericdetector.c squarepixels.c eventfile.c	\
		     impactlistfile.c htrsdetector.c gendet.c		\
		     gendet_init.c genpixgrid.c gendetline.c		\
		     eventlistfile.c ladeventlistfile.c			\
		     genpatternfile.c clocklist.c badpixmap.c		\
		     geneventgrading.c htrseventfile.c			\
		     hexagonalpixels.c arcpixels.c telemetrypacket.c	\
		     htrstelstream.c comadetector.c comaeventfile.c	\
		     psf.c vignetting.c codedmask.c attitudecatalog.c	\
		     attitudefile.c sixt.c vector.c pointsources.c	\
		     extendedsources.c sourceimage.c photon.c		\
		     spectrum.c fits_pha.c check_fov.c			\
		     photonlistfile.c linlightcurve.c rmf.c arf.c	\
		     kdtreeelement.c sourcecatalog.c source.c		\
		     linkedpholist.c erodetbkgrndgen.c phgen.c phimg.c	\
		     phdet.c phproj.c event.c ladevent.c lad.c		\
		     lad_init.c


############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS = genericdetector.h squarepixels.h hexagonalpixels.h	\
		  gaussianchargecloud.h exponentialchargecloud.h	\
		  arcpixels.h eventfile.h impactlistfile.h impact.h	\
		  gendet.h genpixgrid.h gendetline.h clocklist.h	\
		  event.h ladevent.h genpattern.h eventlistfile.h	\
		  ladeventlistfile.h genpatternfile.h badpixmap.h	\
		  geneventgrading.h htrsdetector.h htrseventfile.h	\
		  htrsevent.h telemetrypacket.h htrstelstream.h		\
		  comadetector.h comaeventfile.h comaevent.h psf.h	\
		  vignetting.h codedmask.h attitudecatalog.h		\
		  attitudefile.h telescope.h sixt.h point.h vector.h	\
		  pointsources.h extendedsources.h sourceimage.h	\
		  photon.h spectrum.h fits_pha.h check_fov.h		\
		  photonlistfile.h linlightcurve.h rmf.h arf.h		\
		  kdtreeelement.h sourcecatalog.h source.h		\
		  linkedpholist.h erodetbkgrndgen.h phgen.h phimg.h	\
		  phdet.h phproj.h lad.h

