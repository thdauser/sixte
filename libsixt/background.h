/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2016 Michael Wille, FAU
*/

#ifndef _BACKGROUND_H_
#define _BACKGROUND_H_ 1

#include <ctype.h>
#include <math.h>
#include <sys/timeb.h>

#include "sixt.h"
#include "simput.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "fitsio.h"

/** Needed for seeding the random number generator with
 * a value that changes sufficiently fast */
extern struct timeb time_struct;

/** Additional structure which allows to set auxiliary values
 *  during the background initialization.
 *
 *  To use this functionality it is mandatory to initialize the
 *  background module via 'bkgInitializeAux'.
 */
typedef struct bkgAux {
  double rate;  // sets the rate of the background event file [1/s].
} bkgAux;

/** Structure which contains all the input data from
 * the background simulation file */
typedef struct backgroundInput {
  fitsfile *inputfptr;
  long numrows;
  double interval;
  double intervalsum;
  gsl_rng *randgen;
  
  double xmin_mm;
  double xmax_mm;
  double ymin_mm;
  double ymax_mm;
  double area_sqcm;

  char *timecolname;
  char *energycolname;
  char *xcolname;
  char *ycolname;

  int timecolnum;
  int energycolnum;
  int xcolnum;
  int ycolnum;
  int numevents;
  
  double eventsperinterval;
  double *hit_xpos;
  double *hit_ypos;
  double *hit_time;
  double *hit_energy;

  size_t *eventlist;

  bkgAux aux;
} backgroundInput;

/** Structure which contains the output data and
 * which will be passed to the calling routine */
typedef struct backgroundOutput {
  int numevents;
  int numhits;
  double *hit_xpos;	    // unit as specified in input file
  double *hit_ypos;
  double *hit_time;
  double *hit_energy;  // is required to be in keV or eV.
} backgroundOutput;

/** Information extracted from a lightcurve which allows to modify
 * the background rate with an arbitrary function */
typedef struct backgroundRateFct {
  long numelements;

  double* time;
  float* rate;

  double starttime;
  double* currenttime;
  double intervalsum;
  float* currentrate;
  float currentslope;
} backgroundRateFct;

struct rateCurrentInterval {
  long numelements;
  double interval;

  float* rate;
  long ratesize;
};

/** Use a Simput light curve to manipulate the background rate over time.
 *  The use of this function is optional.
 */
void bkgSetRateFct(const char* const filename, int* const status);

/** Return a randomly chosen eventlist for the given interval;
 * the number of events in the list is determined by poisson statistics
 * while the choice of the events themselves is based on a flat random
 * distribution. This function requires bkgInitialize to be called
 * once before being ready to deliver data. As soon as no more data will
 * be needed one should call the function bkgCleanUp in order to
 * release memory and close the input files. */
backgroundOutput* bkgGetBackgroundList(const double interval);

/* Initialize data structures and read background data table.
 * This function needs to be called once before background events can
 * be requested.
 */
void bkgInitialize(const char* const filename, const unsigned int seed, int* const status);

/** Wrapper function for bkgInitialize() which allows user specified rates
 *  for the background file.
 */
void bkgInitialize_Rate(const char* const filename,
    const unsigned int seed,
    const unsigned int rate,
    int* const status);

/** Wrapper function for bkgInitialize() which allows user specified rates
 *  for the background file with the bkgAux-structure.
 */
void bkgInitializeAux(const char* const filename,
    const unsigned int seed,
    bkgAux* bkgaux,
    int* const status);

/** Free memory of passed backgroundOutput structure */
void bkgFree(backgroundOutput* struct_to_free);

/** Try to cleanup everything and close the input file */
void bkgCleanUp(int* const status);

/************** DEPRECATED section ****************/

/* DEPRECATED: defined for compatibility reasons. */
DEPRECATED(typedef struct backgroundInput eroBackgroundInput);

/* DEPRECATED: defined for compatibility reasons. */
DEPRECATED(typedef struct backgroundOutput eroBackgroundOutput);

/* DEPRECATED: defined for compatibility reasons. */
DEPRECATED(typedef struct backgroundRateFct eroBackgroundRateFct);

/** DEPRECATED: see bkgSetRateFct. */
DEPRECATED(void eroBkgSetRateFct(const char* const filename, int* const status));

/** DEPRECATED: see bkgGetBackgroundList. */
DEPRECATED(eroBackgroundOutput* eroBkgGetBackgroundList(const double interval));

/** DEPRECATED: see bkgInitialize. */
DEPRECATED(void eroBkgInitialize(const char* const filename, const unsigned int seed, int* const status));

/** DEPRECATED: see bkgFree. */
DEPRECATED(void eroBkgFree(backgroundOutput* struct_to_free));

/** DEPRECATED: see bkgCleanUp. */
DEPRECATED(void eroBkgCleanUp(int* const status));

/************* End DEPRECATED section ***************/

#endif   /* _ERODETBKGRNDGEN_H_ */
