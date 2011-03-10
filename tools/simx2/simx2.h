#ifndef SIMX2_H
#define SIMX2_H 1

#include "sixt.h"

#include <wcslib/wcslib.h>

#define TOOLSUB simx2_main
#include "headas_main.c"


struct Parameters {
  char xml_filename[MAXMSG];

  char attitude_filename[MAXMSG];
  float pointing_ra, pointing_dec;

  char simput_filename[MAXMSG];

  char eventlist_filename[MAXMSG];

  double t0, exposure;

  int random_seed;
  
  char fits_templates[MAXMSG];
};


int simx2_getpar(struct Parameters* const par);


#endif /* SIMX2_H */

