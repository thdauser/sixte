#ifndef PHOTON_GENERATION_H
#define PHOTON_GENERATION_H 1

#include "sixt.h"

#include <gsl/gsl_randist.h>

#include "sourceimage.h"
#include "vector.h"
#include "spectrum.h"
#include "photon.h"
#include "astrosources.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "detectors.h"
#include "check_fov.h"


#define TOOLSUB photon_generation_main
#include "headas_main.c"


#endif /* PHOTON_GENERATION_H */

