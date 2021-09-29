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

#ifndef MXS_H
#define MXS_H 1

#include "simput.h"
#include "sixt.h"
#include "linkedimplist.h"
#include "advdet.h"
#include "point.h"
#include <gsl/gsl_rng.h>

// Source ID associated with mxs photons
#define MXS_SRC_ID -10
#define MXS_PH_ID -10

typedef struct {
  double* energy; // [keV]
  long n_events;
} MXSEventList;

typedef struct {
  double mxs_frequency;        // Frequency of the mxs flashes (Hz)
  double mxs_flash_duration;   // Duration of the mxs flashes (s)
  double mxs_rate_det;         // MXS count rate on detector during flash (cps)

  MXSEventList* mxs_eventlist; // Energies MXS photons for sampling
  SimputImg* mxs_img;          // Image of the spatial distribution of the
                               // observed flux.
} MXSparams;

/**
 * Creates a new MXSparams struct. Loads and sets all members and returns
 * a pointer to the new MXSparams struct.
 */
MXSparams* loadMXSparams(double mxs_frequency, double mxs_flash_duration,
                         double mxs_rate_det, char* mxs_filename,
                         int* const status);

/**
 * Frees a previously allocated MXSparams struct.
 */
void freeMXSParams(MXSparams** mxs_params);

/**
 * Generates the next mxs photon that hits the detector. Adjusts flash_start_time
 * and flash_end_time if necessary.
 */
int phmxsgen(double tend, Impact* impact,
             MXSparams* mxs_params, AdvDet* det, double *flash_start_time,
             double *flash_end_time, int* const status);

/**
 * Given the impact time of the previous mxs photon, calculates the impact time
 * of the next mxs photon (mxs photon impact times are modeled as a Poisson
 * process).
 */
double getNextImpactTime(double prevtime, double mxs_rate_det, int* const status);

#endif /* MXS_H */
