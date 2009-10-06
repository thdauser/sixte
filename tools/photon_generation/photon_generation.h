#ifndef PHOTON_GENERATION_H
#define PHOTON_GENERATION_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include <wcs.h>
#include <wcshdr.h>
//#include <wcsfix.h>

#include "sourceimage.h"
#include "vector.h"
#include "spectrum.h"
#include "photon.h"
#include "photonlistfile.h"
#include "astrosources.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "genericdetector.h"
#include "check_fov.h"


#define TOOLSUB photon_generation_main
#include "headas_main.c"


struct Parameters {
  char attitude_filename[MAXMSG];
  char rmf_filename[MAXMSG];
  char sourcelist_filename[MAXMSG];
  char photonlist_filename[MAXMSG];
  char photonlist_template[MAXMSG];

  SourceCategory source_category;

  double t0, timespan;
  double fov_diameter;
};


#endif /* PHOTON_GENERATION_H */

