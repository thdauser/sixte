#ifndef MXS_H
#define MXS_H 1

#include "sixt.h"
#include "linkedimplist.h"
#include "advdet.h"
#include "point.h"
#include <gsl/gsl_rng.h>

// Source ID associated with mxs photons
#define MXS_SRC_ID -10
#define MXS_PH_ID -10

typedef struct {
  double mxs_frequency;      // Frequency of the mxs flashes (Hz)
  double mxs_flash_duration; // Duration of the mxs flashes (s)
  double mxs_rate_det;       // MXS count rate on detector during flash (cps)
} MXSparams;

/**
 * Creates a new MXSparams struct. Loads and sets all members and returns
 * a pointer to the new MXSparams struct.
 */
MXSparams* loadMXSparams (double mxs_frequency, double mxs_flash_duration,
                          double mxs_rate, int numpix, unsigned int seed,
                          int* const status);

/**
 * Generates the next mxs photon that hits the detector. Adjusts flash_start_time
 * and flash_end_time if necessary.
 */
int phmxsgen(AdvDet *det, double tend, Impact* impact,
             MXSparams *mxs_params, double *flash_start_time,
             double *flash_end_time, int* const status);

/**
 * Given the impact time of the previous mxs photon, calculates the impact time
 * of the next mxs photon (mxs photon impact times are modeled as a Poisson
 * process).
 */
double getNextImpactTime(double prevtime, double mxs_rate_det, int* const status);

/**
 * Returns energy of an MXS photon.
 * For 2/3 of the photons, draw a random energy between 3 and 12 keV
 * (Brehmstrahlung continuum), for 1/6 of the photons, emit a 5 keV photons
 * (~ Cr line) and for the last 1/6, emit an 8 keV photon.
 */
double getMXSEnergy();

#endif /* MXS_H */
