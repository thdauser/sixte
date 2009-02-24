#ifndef GENERATE_PHOTONS_H
#define GENERATE_PHOTONS_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
//#include <limits.h>
#include <gsl/gsl_randist.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "headas_rand.h"

#include "sixt.h"
//#include "vector.h"
#include "spectrum.h"
//#include "photon.h"
#include "astrosources.h"
#include "telescope.h"
#include "orbatt.h"
#include "detectors.h"

#define TOOLSUB generate_photons_main
#include "headas_main.c"



#endif /* GENERATE_PHOTONS_H */

