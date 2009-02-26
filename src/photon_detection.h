#ifndef PHOTON_DETECTION_H
#define PHOTON_DETECTION_H 1


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"
#include "detectors.h"
#include "photon.h"
#include "eventlist.h"
#include "point.h"


#define TOOLSUB photon_detection_main
#include "headas_main.c"


//#define MAX(x,y) ( (x)<=(y) ?(y) :(x) )
//#define MIN(x,y) ( (x)<=(y) ?(x) :(y) )



////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int photon_detection_getpar(char impactlist_filename[], char rmf_filename[], 
			    char eventlist_filename[], 
			    double *t0, double *timespan, struct Detector *);


#endif /* PHOTON_DETECTION_H */

