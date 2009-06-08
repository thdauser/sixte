/** 
 * This file is the header file of 'orbatt.c'
 * contains all definitions for orbit and attitude handling.
 */

#ifndef ATTITUDECATALOG_H
#define ATTITUDECATALOG_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"
#include "attitudefile.h"
#include "vector.h"
#include "telescope.h"



typedef struct {
  /** Point of time for which this attitude is valid. */
  double time;
  /** Telescope pointing direction. */
  Vector nz;
  /** Defines the detector x-direction.
   * The x-axis doesn't necessarily have to point in the direction of the telescope
   * motion, but can be distorted by the roll-angle. */
  Vector nx;
} AttitudeEntry;

typedef struct {
  /** Number of AttituideEntry elements in the AttitudeCatalog. */
  long nentries;
  /** Individual AttitudeEntry elements giving the attitude of the telescope
   * at a particular point of time. */
  AttitudeEntry* entry;
} AttitudeCatalog;



/** Constructor for the AttitudeCatalog. */
AttitudeCatalog* get_AttitudeCatalog(const char attitude_filename[],
				     double t0, double timespan, int* status);

/** Destructor for the AttitudeCatalog. */
void free_AttitudeCatalog(AttitudeCatalog* ac);


#endif /* ATTITUDECATALOG_H */

