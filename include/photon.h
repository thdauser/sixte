/** 
 * Header file of "photon.c".
 * It contains all definitions for photon handling.
 */

#ifndef PHOTON_H
#define PHOTON_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

// GSL header files
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_randist.h>


#include "headas.h"
#include "headas_error.h"

#include "vector.h"
#include "detectors.types.h"
#include "astrosources.types.h"


long photon_counter;

// Constants defining the light curve arrays
#define N_LIGHTCURVE_BINS 2048
#define LIGHTCURVE_BINWIDTH 0.0005

// The light curve is constant, without red noise 
// (so the photons have simply Poisson distribution).
//#define CONSTANT_LIGHTCURVE 1   
                                

// The following macros are used to the store light curve and the PSD 
// in the right format for the GSL routines.
#define REAL(z,i) ((z)[(i)])
#define IMAG(z,i) ((z)[N_LIGHTCURVE_BINS-(i)])


// Structure containing all information about a single photon in the sky
struct Photon {
  struct vector direction; // direction from which the photon originates 
                           // (source direction)
  float energy;            // photon energy in [keV]
  double time;             // real time, when the photon is falling on the detector
};



// Structure containing a photon and a pointer to the next photon in the 
// time-ordered photon list.
struct Photon_Entry {
  struct Photon photon;             // photon
  struct Photon_Entry *next_entry;  // pointer to the next entry
};



// Structure representing a bin in the lightcurve.
struct lightcurve_entry {
  double t;           // lower time of bin
  double rate;        // source count rate within the time bin
};


// Include own header files.
#include "detectors.h"
#include "astrosources.h"
#include "random.h"


//////////////////////////////////////////////////////////////////////////
//   Function declarations
//////////////////////////////////////////////////////////////////////////


// Creates photons according to a particular rate specified by the given light curve
// and adds them to the time ordered photon list.
// The return value is the value of the error status variable.
int create_photons(struct source_cat_entry* src, double time, double dt,
		   struct Photon_Entry** pl, struct Detector detector, 
		   gsl_rng *gsl_random_g);

// Clears the photon list.
void clear_photon_list(struct Photon_Entry **);

// Inserts a new photon into the time ordered photon list.
int insert_photon(struct Photon_Entry **, struct Photon);

// Creates a randomly chosen photon energy according to the spectrum of the 
// specified source.
float photon_energy(struct source_cat_entry src, struct Detector);

// Function produces a light curve for a given source.
int create_lightcurve(struct source_cat_entry *src, double time, gsl_rng *gsl_random_g);


// This function creates a new event list table in the specified FITS file.
// It also inserts  header information.
// The function returns '0', if it is run successfully.
// Otherwise the return value is '1'.
int create_photonlist_file(struct Eventlist_File*, struct Detector detector,
			   double tstart, double tend, char *telescope_name,
			   char *ccd_name, char *instrument_name,
			   int *status);


// Open an existing FITS file and try to get the first extension that contains
// a binary table.
int open_photonlist_file(struct Eventlist_File*, int* status);


#endif

