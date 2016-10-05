#ifndef TESGENPIXIMP_H
#define TESGENPIXIMP_H

#include "sixt.h"

#include "sixteconfig.h"
#include "pixelimpactfile.h"

#include <gsl/gsl_rng.h>

#include <assert.h>
#include <parinput.h>
#include <stdlib.h>
#include <math.h>

// generation modes
#define MODE_CONST 1
#define MODE_LIN 2
#define MODE_RAND 3
#define MODE_SINRAND 4


// input parameter struct
typedef struct {
  char *pixFilename; // filename of piximplist
  int mode;          // generation mode
  int pixid;         // pixel ID
  double tstart;     // start time
  double tstop;      // stop time
  double dtau;       // time delay between photons
  double EConst;     // constant photon energy
  double Emin;       // minimum photon energy
  double Emax;       // maximum photon energy
  int nPhotons;   // nr of photons to generate
  long seed;     // seed for RNG
  double offset;     // constant offset of sin distribution
  double amplitude;  // amplitude of sin distribution
  double f;          // frequency of sin distribution
  double shift;      // phase-shift of sin distribution
} tesgenimppars;

// parameter input function
void tesgenimpacts_getpar(tesgenimppars *par, int *status);

// comparision function, for qsort
int cmp(const void *x, const void *y);

// sine probability distribution. Not normalized!
double sin_dist(double t, double offset, double amplitude, double f, double shift);

// indefinite integral of sin_dist
double sin_int(double t, double offset, double amplitude, double f, double shift);

// get a photon impact time from the sine distribution using
// an acceptance-rejection method
double sin_get_impact_time(gsl_rng *rng, double offset, double amplitude, double f, double tstart, double tstop, double shift);

#endif /* TESGENPIXIMP_H */
