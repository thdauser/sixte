#include "mxs.h"

// Random number generator
gsl_rng *rng = NULL; // initialize to NULL and set it in loadMXSparams

/**
 * Creates a new MXSparams struct. Loads and sets all members and returns
 * a pointer to the new MXSparams struct.
 */
MXSparams* loadMXSparams (double mxs_frequency, double mxs_flash_duration,
                          double mxs_rate, int numpix, unsigned int seed,
                          int* const status)
{
  MXSparams *mxs_params = malloc(sizeof(MXSparams));
  CHECK_NULL(mxs_params, *status, "Memory allocation for mxs_params failed");

  mxs_params->mxs_frequency = mxs_frequency;
  mxs_params->mxs_flash_duration = mxs_flash_duration;
  // Multiply mxs_rate by numpix to get the mxs rate on the whole detector.
  mxs_params->mxs_rate_det = mxs_rate * numpix;

  // Initialize the rng
  if (rng == NULL) {
    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, seed);
  }

  return mxs_params;
}

/**
 * Generates the next mxs photon that hits the detector. Adjusts flash_start_time
 * and flash_end_time if necessary.
 */
int phmxsgen(AdvDet *det, double tend, Impact* impact,
             MXSparams *mxs_params, double *flash_start_time,
             double *flash_end_time, int* const status)
{
  // Keeps track internally of the last mxs impact time. Used to calculate the
  // next impact time.
  static double time_of_last_mxs_impact = -1.;

  // If phmxsgen is called for the first time in a gti,
  // time_of_last_mxs_impact might be smaller than the start time of the first
  // flash in that gti. Set time_of_last_mxs_impact to the flash_start_time in
  // that case.
  if (time_of_last_mxs_impact < *flash_start_time) {
    time_of_last_mxs_impact = *flash_start_time;
  }

  // Calculate time of the next mxs impact on the detector.
  double time_of_next_mxs_impact = getNextImpactTime(time_of_last_mxs_impact,
                                                     mxs_params->mxs_rate_det,
                                                     status);

  // Check if time_of_next_mxs_impact is still within the current flash. If it is
  // not, set time_of_next_mxs_impact to the first impact during the next flash and
  // calculate new flash interval.
  // (while loop because we could have a flash interval without an mxs photon if
  // getNextImpactTime provides an impact time outside this interval)
  while (time_of_next_mxs_impact > *flash_end_time) {

    time_of_last_mxs_impact = *flash_end_time + 1./mxs_params->mxs_frequency;
    *flash_start_time = time_of_last_mxs_impact;
    *flash_end_time = *flash_start_time + mxs_params->mxs_flash_duration;
    time_of_next_mxs_impact = getNextImpactTime(time_of_last_mxs_impact,
                                                mxs_params->mxs_rate_det,
                                                status);
  }

  if (time_of_next_mxs_impact <= tend) {
    // Get energy of the mxs impact.
    double energy = getMXSEnergy();

    // Generate a random impact location on the detector.
    // First get a random pixel of the array.
    unsigned long int pixid = gsl_rng_uniform_int(rng, det->npix);

    // Get corresponding values of x and y.
    double x_coord = det->sx + det->pix[pixid].sx
                   + (sixt_get_random_number(status)-0.5)*det->pix[pixid].width;

    double y_coord = det->sy + det->pix[pixid].sy
                   + (sixt_get_random_number(status)-0.5)*det->pix[pixid].height;


    // Assign calculated values to the impact
    impact->time = time_of_next_mxs_impact;
    impact->energy = energy;
    impact->position.x = x_coord;
    impact->position.y = y_coord;
    impact->ph_id = MXS_PH_ID;
    impact->src_id = MXS_SRC_ID;

    // Update time_of_last_mxs_impact.
    time_of_last_mxs_impact = time_of_next_mxs_impact;

    // Yes, we have generated an mxs photon.
    return 1;
  } else {
    // Next mxs photon would be outside current gti.
    return 0;
  }
}

/**
 * Given the impact time of the previous mxs photon, calculates the impact time
 * of the next mxs photon (mxs photon impact times are modeled as a Poisson
 * process).
 */
double getNextImpactTime(double prevtime, double mxs_rate_det, int* const status)
{
  // Time intervals between subsequent mxs photons are exponentially distributed.
  return prevtime + (-log(1-sixt_get_random_number(status))/mxs_rate_det);
}

/**
 * Returns energy of an MXS photon.
 * For 2/3 of the photons, draw a random energy between 3 and 12 keV
 * (Brehmstrahlung continuum), for 1/6 of the photons, emit a 5 keV photons
 * (~ Cr line) and for the last 1/6, emit an 8 keV photon.
 */
double getMXSEnergy()
{
  unsigned long int random_int = gsl_rng_uniform_int(rng, 6);

  if (random_int <= 3) {
    return (12 - 3) * gsl_rng_uniform(rng) + 3;
  } else if ( random_int==4 ) {
    return 5.0;
  } else {
    return 8.0;
  }
}
