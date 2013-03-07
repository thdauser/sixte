AM_CFLAGS =-I@top_srcdir@/extlib/simput/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/simput/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp

if RCL
AM_CFLAGS+=-I$(RCL_INC) -DUSE_RCL
AM_LDFLAGS=-L$(RCL_LIB) -lrcl
endif

########## DIRECTORIES ###############

# Directory where to install the FITS template files.
fitstemplatesdir = $(pkgdatadir)/templates
dist_fitstemplates_DATA =coma.eventlist.tpl
dist_fitstemplates_DATA+=eventlist.tpl
dist_fitstemplates_DATA+=ladsignallist.tpl
dist_fitstemplates_DATA+=ladeventlist.tpl
dist_fitstemplates_DATA+=ladspec.tpl
dist_fitstemplates_DATA+=genlc.tpl
dist_fitstemplates_DATA+=eroeventlist.tpl
dist_fitstemplates_DATA+=patternlist.tpl
dist_fitstemplates_DATA+=htrs.eventlist.tpl
dist_fitstemplates_DATA+=impactlist.tpl
dist_fitstemplates_DATA+=ladimpactlist.tpl
dist_fitstemplates_DATA+=photonlist.tpl 
dist_fitstemplates_DATA+=gti.tpl

# Directory where to install the XML files.
xmldir=$(pkgdatadir)/xml
dist_xml_DATA =wfi.iaat.xml
dist_xml_DATA+=htrs.xml
dist_xml_DATA+=comadet.xml
dist_xml_DATA+=aciss3.xml

AM_CFLAGS+=-DSIXT_DATA_PATH='"$(pkgdatadir)"'

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libsixt.la

# Sources:
libsixt_la_SOURCES=genericdetector.c squarepixels.c eventfile.c		\
		   impactlistfile.c ladimpactlistfile.c			\
		   htrsdetector.c geninst.c gendet.c gentel.c		\
		   genpixgrid.c gendetline.c eventlistfile.c		\
		   ladsignallistfile.c ladeventlistfile.c phabkg.c	\
		   patternfile.c clocklist.c badpixmap.c		\
		   htrseventfile.c hexagonalpixels.c arcpixels.c	\
		   telemetrypacket.c htrstelstream.c comadetector.c	\
		   comaeventfile.c psf.c vignetting.c codedmask.c	\
		   attitudecatalog.c attitudefile.c sixt.c vector.c	\
		   photon.c fits_pha.c check_fov.c photonlistfile.c	\
		   rmf.c arf.c kdtreeelement.c sourcecatalog.c		\
		   source.c linkedpholist.c linkedladsiglist.c		\
		   erodetbkgrndgen.c phgen.c phimg.c phdet.c phproj.c	\
		   phpat.c event.c pattern.c ladsignal.c ladevent.c	\
		   ladimpact.c lad.c lad_init.c xmlbuffer.c		\
		   spectrum.c gti.c sourceimage.c det_phi_max.c		\
		   attitudecatalog2.c reconstruction.c

libsixt_la_LIBADD =@top_builddir@/extlib/heainit/libhdinit.la
libsixt_la_LIBADD+=@top_builddir@/extlib/heaio/libhdio.la
libsixt_la_LIBADD+=@top_builddir@/extlib/heautils/libhdutils.la
libsixt_la_LIBADD+=@top_builddir@/extlib/heasp/libhdsp.la
libsixt_la_LIBADD+=@top_builddir@/extlib/ape/src/libape.la
#libsixt_la_LIBADD+=@top_builddir@/extlib/simput/libsimput/libsimput.la
libsixt_la_LIBADD+=@top_builddir@/extlib/simput/extlib/cfitsio/libcfitsio.a

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=genericdetector.h squarepixels.h hexagonalpixels.h	\
		gaussianchargecloud.h exponentialchargecloud.h		\
		arcpixels.h eventfile.h impactlistfile.h		\
		ladimpactlistfile.h impact.h geninst.h gendet.h		\
		gentel.h genpixgrid.h gendetline.h clocklist.h		\
		event.h ladsignal.h ladevent.h ladimpact.h pattern.h	\
		eventlistfile.h ladsignallistfile.h phabkg.h		\
		ladeventlistfile.h patternfile.h badpixmap.h		\
		htrsdetector.h htrseventfile.h htrsevent.h		\
		telemetrypacket.h htrstelstream.h comadetector.h	\
		comaeventfile.h comaevent.h psf.h vignetting.h		\
		codedmask.h attitudecatalog.h attitudefile.h		\
		telescope.h sixt.h point.h vector.h photon.h		\
		fits_pha.h check_fov.h photonlistfile.h rmf.h arf.h	\
		kdtreeelement.h sourcecatalog.h source.h		\
		linkedpholist.h linkedladsiglist.h erodetbkgrndgen.h	\
		phgen.h phimg.h phdet.h phproj.h phpat.h lad.h		\
		xmlbuffer.h spectrum.h gti.h sourceimage.h		\
		det_phi_max.h attitudecatalog2.c reconstruction.h			

