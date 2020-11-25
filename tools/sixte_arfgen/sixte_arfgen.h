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

   Copyright 2020 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#ifndef SIXTE_ARFGEN_H
#define SIXTE_ARFGEN_H 1

#include "sixt.h"
#include "parinput.h"
#include "geninst.h"
#include "photon.h"
#include "attitude.h"
#include "phimg.h"
#include "phproj.h"
#include "eventfile.h"
#include "radec2xylib.h"
#include <gsl/gsl_interp.h>
#include <stdbool.h>

const int sampling_factor = 10; // Sample every 10-th ARF bin

/**
 * Tool parameters
 */
typedef struct {
  char* ARFCorr;         // Corrected ARF output file
  char* XMLFile;         // XML input file with instrument definition
  int Seed;              // Seed for RNG
  int clobber;           // Overwrite output files if exist?

  double PointingRA;     // Right ascension of telescope pointing [degree]
  double PointingDec;    // Declination of telescope pointing [degree]
  double rollangle;      // Roll angle of telescope pointing [degree]
  char* Attitude;        // Attitude input file
  double SourceRA;       // Right ascension of source [degree]
  double SourceDec;      // Declination of source [degree]
  char* ImageFile;       // Count image (containing WCS definition) used for region definition
  char Projection[MAXFILENAME]; // WCS projection type (usually SIN)
  float RefRA;           // Right ascension of WCS reference point [deg]
  float RefDec;          // Declination of WCS reference point [deg]
  char regfilter[MAXFILENAME]; // Region filter

  double MJDREF;         // Reference Modified Julian Date
  double TSTART;
  double exposure;

  int n_photons;         // Number of photons simulated per ARF bin
  double photon_rate;    // Photon rate [cts/s]

  char eventlist_filename[MAXFILENAME]; // Filename of temporary eventlist
} Parameters;


typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;

  /** Column numbers. */
  int cenerg_lo, cenerg_hi, cspecresp;
} ARFFile;

typedef struct {
  int* bin_idxs;       // ARF bins for which correction factors are calculated
  double* corr_fac;    // Associated correction factors
  double* corr_energ;  // for specified energies (bin midpoints).
  int n_corr;          // Total number of correction factors
} ARFCorr;

////////////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////////////

// Initializes parameters struct and loads all data members.
int sixte_arfgen_getpar(Parameters* const par);

// Frees parameters struct and all its data members.
void sixte_arfgen_freepar(Parameters* par);

// Generates the corrected ARF.
// For every n-th ARF bin (where n is given by sampling_factor), n_photons photons
// (originating from extraction region) are generated and then imaged.
// ARF correction factors are calculated for these bins, given by
// n_impacts/n_photons. Multiplication with the original effective area of the
// corresponding bins gives the corrected effective areas. The whole ARF is then
// calculated by interpolation.
void sixte_arfgen(GenInst* const inst, EventFile* ilf,
            Attitude* const ac, ARFFile* arf_out,
            Parameters* par, int* status);

// Calculates the correction factors (n_impacts/n_photons) for every n-th ARF
// bin (where n is given by sampling_factor).
void calc_arf_corr(ARFCorr* arf_corr, const struct ARF* const arf,
                   GenInst* const inst, Attitude* const ac, EventFile* elf,
                   Parameters* par, int* status);

// Given an array of correction factors (for every n-th ARF bin) and corresponding
// energies, updates effective areas in arf_out by interpolation.
void update_arf(ARFFile* arf_out, ARFCorr* arf_corr,
                const struct ARF* const arf, int* status);

// Generates new ARF file, initialized as a copy of arf_filepathname.
ARFFile* initARFFile(const Parameters* const par, char* arf_filepathname, int* status);

// Initialize the ARF correction factors. Currently, every sampling_factor-th
// bin is simulated (last bin is always included to avoid extrapolation when writing
// the corrected ARF).
ARFCorr* initARFCorr(const struct ARF* const arf, int* status);

// Frees all memory related to arf_corr
void freeARFCorr(ARFCorr** arf_corr);

// Frees an ARF file
void freeARFFile(ARFFile** arf_file, int* status);

// Generates photons with energies in [E_low, E_high) at (RA, Dec) and given rate.
void gen_photon(Photon* ph, float E_low, float E_high, double RA,
                double Dec, double time, int* status);

// Sets the energy of a photon to random value in [E_low, E_high)
void set_ph_energy(Photon* ph, double E_low, double E_high, int* status);

// Sets the origin of a photon to (RA, Dec).
void set_ph_position(Photon* ph, double RA, double Dec);

// Sets photon meta info: time (considering given rate), id (ascending)
// and src id (=0).
void set_ph_info(Photon* ph, double rate);

#endif /* SIXTE_ARFGEN_H */
