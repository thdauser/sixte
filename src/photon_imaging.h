#ifndef IMAGING_H
#define IMAGING_H 1


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>
//#include <gsl/gsl_randist.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"
#include "vector.h"
//#include "spectrum.h"
#include "check_fov.h"
#include "fits_ctlg.h"
//#include "imglib.h"
#include "random.h"
//#include "strftcpy.h"
#include "psf.h"
#include "photon.h"
//#include "sources.h"
#include "telescope.h"
#include "orbatt.h"


#define TOOLSUB photon_imaging_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int photon_imaging_getpar(char photonlist_filename[], char orbit_filename[], 
			  char attitude_filename[], 
			  char psf_filename[], double *t0, double *timespan, 
			  struct Telescope *, struct Detector *);


#endif
