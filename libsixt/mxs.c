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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include "mxs.h"


static MXSEventList* loadMXSEvents(char* mxs_filename, int* const status) {
  MXSEventList* mxs_eventlist = malloc(sizeof(MXSEventList));
  CHECK_MALLOC_RET_NULL_STATUS(mxs_eventlist, *status)

  // Add extension name to filename
  char mxs_filename_ext[MAXFILENAME];
  strcpy(mxs_filename_ext, mxs_filename);
  strcat(mxs_filename_ext, "[EVENTS]");

  // Open MXS file and load all data from EVENTS extension
  fitsfile *fptr;
  fits_open_table(&fptr, mxs_filename, READONLY, status);

  // Get n_events
  fits_get_num_rows(fptr, &(mxs_eventlist->n_events), status);
  CHECK_STATUS_RET(*status, NULL);

  int colnum, anynul;
  // Initialize energy
  mxs_eventlist->energy = malloc(mxs_eventlist->n_events*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(mxs_eventlist->energy, *status);
  fits_get_colnum(fptr, CASEINSEN, "ENERGY", &colnum, status);
  fits_read_col(fptr, TDOUBLE, colnum, 1, 1, mxs_eventlist->n_events, NULL,
                mxs_eventlist->energy, &anynul, status);

  fits_close_file(fptr, status);

  return mxs_eventlist;
}


static SimputImg* loadMXSImage(char* mxs_filename, int* const status) {
  // Add extension name to filename
  char mxs_filename_ext[MAXFILENAME];
  strcpy(mxs_filename_ext, mxs_filename);
  strcat(mxs_filename_ext, "[ILLUMINATION_MAP]");

  return loadSimputImg(mxs_filename_ext, status);
}


/**
 * Creates a new MXSparams struct. Loads and sets all members and returns
 * a pointer to the new MXSparams struct.
 */
MXSparams* loadMXSparams(double mxs_frequency, double mxs_flash_duration,
                         double mxs_rate_det, char* mxs_filename,
                         int* const status) {
  MXSparams *mxs_params = malloc(sizeof(MXSparams));
  CHECK_NULL(mxs_params, *status, "Memory allocation for mxs_params failed");

  mxs_params->mxs_frequency = mxs_frequency;
  mxs_params->mxs_flash_duration = mxs_flash_duration;
  mxs_params->mxs_rate_det = mxs_rate_det;

  mxs_params->mxs_eventlist = loadMXSEvents(mxs_filename, status);
  mxs_params->mxs_img = loadMXSImage(mxs_filename, status);

  return mxs_params;
}


void freeMXSParams(MXSparams** mxs_params) {
  if (*mxs_params != NULL) {
    if ((*mxs_params)->mxs_eventlist != NULL) {
      if ((*mxs_params)->mxs_eventlist->energy != NULL) {
        free((*mxs_params)->mxs_eventlist->energy);
      }

      free((*mxs_params)->mxs_eventlist);
    }

    if ((*mxs_params)->mxs_img != NULL) {
      free((*mxs_params)->mxs_img);
    }

    free(*mxs_params);
    *mxs_params = NULL;
  }
}


/**
 * Given a sorted array of Edeps, samples an energy deposition of a cosmic from
 * this array, using the (empirical) cumulative distribution function and
 * inverse transform sampling (see Nelson B., 2013, Foundations and Methods of
 * Stochastic Simulation: A First Course, p. 118).
 */
static float getMXSEnergy(double *energies, long n_events, int* const status) {
  // Generate random number in (0,1)
  double u;
  do {
    u = sixt_get_random_number(status);
  } while (u == 0);

  long i = (long)ceil( (n_events-1)*u );
  return (float) (energies[i-1] + (n_events-1) * (energies[i]-energies[i-1])
                                               * (u - (i-1.)/(n_events-1.)));
}


static int positionOnPixelArray(AdvDet* det, struct Point2d* position) {
  // Transform impact coordinates into detector coordinate system
  Impact detimp;
  detimp.position.x = position->x - det->sx;
  detimp.position.y = position->y - det->sy;

  // Check for hit
  for (int ii = 0; ii < det->npix; ii++) {
    if ( CheckAdvPixImpact(det->pix[ii], &detimp) != 0) {
      return 1;
    }
  }

  // Impact not on pixel array
  return 0;
}


static struct Point2d getMXSPosition(AdvDet* det, SimputImg* mxs_img,
                                     int* const status) {
  struct Point2d mxs_position;

  // Sample a random position on the pixel array
  do {
    // Get a random position from the image [pixel].
    double xd, yd;
    drawRndPosFromImg(mxs_img, &xd, &yd, status);

    // Get corresponding impact position [m]
    mxs_position.x = (xd - mxs_img->wcs->crpix[0]) * mxs_img->wcs->cdelt[0];
    mxs_position.y = (yd - mxs_img->wcs->crpix[1]) * mxs_img->wcs->cdelt[1];
  } while ( positionOnPixelArray(det, &mxs_position) == 0 );

  return mxs_position;
}


static Impact genMXSImpact(MXSEventList* mxs_eventlist,
                           SimputImg* mxs_img,
                           AdvDet* det,
                           int* const status) {
  Impact mxs_impact;

  mxs_impact.energy = getMXSEnergy(mxs_eventlist->energy, mxs_eventlist->n_events,
                                   status);
  mxs_impact.position = getMXSPosition(det, mxs_img, status);

  mxs_impact.ph_id = MXS_PH_ID;
  mxs_impact.src_id = MXS_SRC_ID;

  return mxs_impact;
}

/**
 * Generates the next mxs photon that hits the detector. Adjusts flash_start_time
 * and flash_end_time if necessary.
 */
int phmxsgen(double tend, Impact* impact,
             MXSparams* mxs_params, AdvDet* det, double *flash_start_time,
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
    *flash_start_time += 1./mxs_params->mxs_frequency;
    *flash_end_time = *flash_start_time + mxs_params->mxs_flash_duration;

    time_of_last_mxs_impact = *flash_start_time;
    time_of_next_mxs_impact = getNextImpactTime(time_of_last_mxs_impact,
                                                mxs_params->mxs_rate_det,
                                                status);
  }

  if (time_of_next_mxs_impact <= tend) {
    Impact mxs_impact = genMXSImpact(mxs_params->mxs_eventlist,
                                     mxs_params->mxs_img,
                                     det, status);
    copyImpact(impact, &mxs_impact);

    // Assign calculated time to the impact
    impact->time = time_of_next_mxs_impact;

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
