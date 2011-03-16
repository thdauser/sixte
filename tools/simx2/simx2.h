#ifndef SIMX2_H
#define SIMX2_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "gendet.h"
#include "impactlistfile.h"
#include "phgen.h"
#include "phimg.h"
#include "photonlistfile.h"
#include "vector.h"
#include "xraysourcecatalog.h"

#define TOOLSUB simx2_main
#include "headas_main.c"


struct Parameters {
  char xml_filename[MAXFILENAME];

  char attitude_filename[MAXFILENAME];
  /** [deg] */
  float pointing_ra, pointing_dec;

  char simput_filename[MAXFILENAME];
  /** Source type can be point (0), flat (1), or image(2) */
  int src_type;
  /** [erg/s/cm**2] */
  float src_flux;
  /** Source position [deg]. Only applicable for point sources. */
  float src_ra, src_dec;
  char src_spectrum_filename[MAXFILENAME];
  /** Energy of a mono energetic source [keV]. */
  float src_mono_energy;
  char src_image_filename[MAXFILENAME];

  char eventlist_filename[MAXFILENAME];

  double t0, exposure;

  int random_seed;
  
  char fits_templates[MAXFILENAME];
};


int simx2_getpar(struct Parameters* const par);


#endif /* SIMX2_H */

