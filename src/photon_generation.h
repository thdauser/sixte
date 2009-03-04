#ifndef PHOTON_GENERATION_H
#define PHOTON_GENERATION_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <limits.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "headas_rand.h"

#include "sixt.h"
#include "vector.h"
#include "spectrum.h"
#include "photon.h"
#include "astrosources.h"
#include "telescope.h"
#include "orbatt.h"
#include "detectors.h"
#include "check_fov.h"


#define TOOLSUB photon_generation_main
#include "headas_main.c"



#endif /* PHOTON_GENERATION_H */

