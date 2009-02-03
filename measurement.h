#ifndef MEASUREMENT_H
#define MEASUREMENT_H 1


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>
#include <png.h>
#include <gsl/gsl_randist.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "vector.h"
#include "spectrum.h"
#include "check_fov.h"
#include "fits_ctlg.h"
#include "imglib.h"
#include "random.h"
#include "strftcpy.h"
#include "detector.h"
#include "psf.h"
#include "photon.h"
#include "sources.h"
#include "global_constants.h"
#include "telescope.h"
#include "orbatt.h"
#include "measurement_array.h"
#include "event_list.h"


#define TOOLSUB measurement_main
#include "headas_main.c"

#define LINELENGTH 1024        // maximum linelength in the ASCII file
#define MAXMSG 256             // maximum length of an output/error message

#define MAX(x,y) ( (x)<=(y) ?(y) :(x) )
#define MIN(x,y) ( (x)<=(y) ?(x) :(y) )



////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int measurement_getpar(char orbit_filename[], char attitude_filename[], 
		       int *n_sourcefiles, 
		       char source_filename[MAX_NSOURCEFILES][FILENAME_LENGTH], 
		       char psf_filename[], char detrsp_filename[], 
		       char spectrum_filename[], char eventlist_filename[], 
		       double *t0, 
		       double *timespan, struct Telescope *, struct Detector *,
		       double *bandwidth, float *background_countrate);

// Function performs a measurement for a specific source inside the FOV 
// with given telescope data (direction of axis, direction of motion).
int measure(struct Photon, struct Detector, struct Telescope, struct PSF_Store,
	    struct Event_List_File* event_list_file);



#endif
