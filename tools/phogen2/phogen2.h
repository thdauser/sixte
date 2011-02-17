#ifndef PHOGEN2_H
#define PHOGEN2_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include <wcslib/wcslib.h>

#include "interface.h"
#include "sixt_string.h"
#include "sourceimage.h"
#include "pointsources.h"
#include "pointsourcefile.h"
#include "pointsourcecatalog.h"
#include "pointsourcelist.h"
#include "extendedsources.h"
#include "vector.h"
#include "spectrum.h"
#include "photon.h"
#include "photonlistfile.h"
#include "astrosources.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "arf.h"
#include "check_fov.h"
#include "kdtree.h"
#include "gendet.h"

#define TOOLSUB phogen_main
#include "headas_main.c"


struct Parameters {
  char xml_filename[MAXMSG];
  char attitude_filename[MAXMSG];
  char simput_filename[MAXMSG];
  char photonlist_filename[MAXMSG];
  char photonlist_template[MAXMSG];

  double t0, timespan;
};


int phogen_getpar(struct Parameters* parameters);


#endif /* PHOGEN2_H */

