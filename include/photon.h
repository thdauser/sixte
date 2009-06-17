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

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "vector.h"
#include "detectors.types.h"


long photon_counter;

// Constants defining the light curve arrays
#define N_LIGHTCURVE_BINS     (2048)
#define LIGHTCURVE_BINWIDTH   (0.0005)
#define N_PHOTON_FIELDS       (4)     // TIME, ENERGY, RA, DEC
#define N_IMPACT_FIELDS       (4)     // TIME, ENERGY, X, Y


// The light curve is constant, without red noise 
// (so the photons have simply Poisson distribution).
#define CONSTANT_LIGHTCURVE 1   
                                

// The following macros are used to the store light curve and the PSD 
// in the right format for the GSL routines.
#define REAL(z,i) ((z)[(i)])
#define IMAG(z,i) ((z)[N_LIGHTCURVE_BINS-(i)])



/** Contains all information about a single photon in the sky. */
typedef struct {
  double time;  /**< Real time, when the photon is falling on the detector (in [s]). */
  float energy; /**< Photon energy in [keV]. */

  double ra, dec; /**< Right ascension and declination of photon position [rad]. */

  // REMOVE
  Vector direction; // direction from which the photon originates 
                    // (source direction)
} Photon;



// Structure containing a photon and a pointer to the next photon in the 
// time-ordered photon list.
struct PhotonOrderedListEntry {
  Photon photon; 
  struct PhotonOrderedListEntry *next;  // pointer to the next entry
};


/** Entry in the binary tree that stores the generated photons. */
struct PhotonBinaryTreeEntry {
  Photon photon; /**< Photon data. */

  struct PhotonBinaryTreeEntry* sptr; /**< Pointer to entry with smaller time value. */
  struct PhotonBinaryTreeEntry* gptr; /**< Pointer to entry with greater time value. */
};


// Structure representing a bin in the lightcurve.
struct lightcurve_entry {
  double t;           // lower time of bin
  double rate;        // source count rate within the time bin
};




// Include own header files.
#include "detectors.h"
#include "pointsources.h"
#include "random.h"


//////////////////////////////////////////////////////////////////////////
//   Function declarations
//////////////////////////////////////////////////////////////////////////


// Creates photons according to a particular rate specified by the given 
// light curve and adds them to the time ordered photon list.
// The return value is the value of the error status variable.
int create_photons(PointSource* ps, double time, double dt,
		   struct PhotonOrderedListEntry** pl, Detector*, gsl_rng *gsl_random_g);

// Clears the photon list.
void clear_photon_list(struct PhotonOrderedListEntry **);

// Inserts a new photon into the time ordered photon list.
int insert_photon(struct PhotonOrderedListEntry **, Photon);

// Creates a randomly chosen photon energy according to the spectrum of the 
// specified source.
//float photon_energy(struct source_cat_entry src, Detector*);
float photon_energy(struct Spectrum*, Detector*);

// Function produces a light curve for a given source.
//int create_lightcurve(struct source_cat_entry *src, double time, 
int create_lightcurve(PointSource* ps, double time, gsl_rng *gsl_random_g);


// This routine creates a new FITS file with a binary table to store a photon list
// of photons from astronomical x-ray sources. The photon list can be read by a 
// telescope simulation for further processing and calculating the impact positions
// of the photons on the detector.
int create_photonlist_file(fitsfile **, char filename[], int *status);


// This routine creates a new FITS file with a binary table to store an impact list
// of photons on the detector. The list can be further processed by a detector
// simulation to create an event list.
int create_impactlist_file(fitsfile **, char filename[], int *status);


/** Insert a new photon to an existing binary tree.
 * The pointer to the PhotonBinaryTreeEntry can be NULL. */
void insert_Photon2BinaryTree(struct PhotonBinaryTreeEntry**, Photon*);


#endif  /* PHOTON_H */

