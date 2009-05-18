/** 
 * This file is the header file of 'orbatt.c'
 * contains all definitions for orbit and attitude handling.
 */

#ifndef ORBATT_H
#define ORBATT_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "fits_attitude.h"
#include "fits_ctlg.h"
#include "vector.h"
#include "telescope.h"


// maximum number of entries in the orbit-/attitude-catalog (satellite_catalog)
//#define MAX_NORBIT_ENTRIES 500000   


typedef struct {
  /** Point of time for which this attitude is valid. */
  double time;
  /** Telescope pointing direction. */
  struct vector nz;
  /** Defines the detector x-direction.
   * The x-axis doesn't necessarily have to point in the direction of the telescope
   * motion, but can be distorted by the roll-angle. */
  struct vector nx;
} AttitudeEntry;

typedef struct {
  /** Number of AttituideEntry elements in the AttitudeCatalog. */
  long nentries;
  /** Individual AttitudeEntry elements giving the attitude of the telescope
   * at a particular point of time. */
  AttitudeEntry* entry;
} AttitudeCatalog;


// This function allocates memory for the satellite catalog 
// (containing orbit and attitude information), reads this data 
// from the corresponding FITS files and stores it in the catalog.
/*int get_satellite_catalog(struct Telescope **, long *nentries, double t0, 
			  double timespan, const char orbit_filename[], 
			  const char attitude_filename[]);*/


/** Constructor for the AttitudeCatalog. */
AttitudeCatalog* get_AttitudeCatalog(const char attitude_filename[],
				     double t0, double timespan, int* status);

/** Destructor for the AttitudeCatalog */
void free_AttitudeCatalog(AttitudeCatalog* ac);


#endif /* ORBATT_H */

