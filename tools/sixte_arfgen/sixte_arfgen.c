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

#include "sixte_arfgen.h"

#define TOOLSUB sixte_arfgen_main
#include "sixt_main.c"

int sixte_arfgen_main()
{
  // Register HEATOOL
  set_toolname("sixte_arfgen");
  set_toolversion("0.01");

  // Program parameters
  Parameters par;

  // Instrument setup
  GenInst* inst = NULL;

  // Attitude
  Attitude* ac = NULL;

  // Temporary event file
  EventFile* elf = NULL;
  EventFile* memelf = NULL;

  // Output ARF
  ARFFile* arf_out = NULL;

  // Error status.
  int status = EXIT_SUCCESS;

  do { // Beginning of ERROR HANDLING Loop

    // ---- Initialization ----

    // Read the parameters using PIL
    status = sixte_arfgen_getpar(&par);

    // Initialize the random number generator
    unsigned int seed = getSeed(par.Seed);
    sixt_init_rng(seed, &status);

    // Load the instrument configuration
    inst = loadGenInst(par.XMLFile, seed, &status);

    // Set up the Attitude
    if (par.Attitude == NULL) {
      // Set up a pointing Attitude
      par.TSTART = 0;
      par.exposure = par.n_photons*(1./par.photon_rate); // Exposure for one bin
      ac = getPointingAttitude(par.MJDREF, par.TSTART, par.TSTART + par.exposure,
                               par.PointingRA*M_PI/180., par.PointingDec*M_PI/180.,
                               par.rollangle*M_PI/180., &status);
      CHECK_STATUS_BREAK(status);
    } else {
      // Load the attitude from the given file.
      ac = loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Adjust the photon rate
      if (par.TSTART >= 0 && par.exposure > 0) { // user requests specific subset of attitude
        // Check if the required time interval is a subset of the period described by the attitude file.
        checkAttitudeTimeCoverage(ac, par.MJDREF, par.TSTART,
          par.TSTART + par.exposure, &status);
        CHECK_STATUS_BREAK(status);

        // Sample photon times uniformly during exposure
        par.photon_rate = par.n_photons / par.exposure;
      } else {
        // Sample complete attitude
        par.photon_rate = par.n_photons / (ac->tstop - ac->tstart);
      }
    }

    // Open temporary event file
    elf = openNewEventFile(par.eventlist_filename, "", "", "",
			                     inst->tel->arf_filename, inst->det->rmf_filename,
		                    	 par.MJDREF, 0.0, par.TSTART, par.TSTART + par.exposure,
		                    	 inst->det->pixgrid->xwidth,
	                    		 inst->det->pixgrid->ywidth,
	                    		 1, &status);
    memelf = copyEventFileMemory(elf, &status);

    // Add EVTYPE keyword (required by radec2xy)
    fits_update_key(memelf->fptr, TSTRING, "EVTYPE", "PATTERN",
                    "event type", &status);

    // Initialize output ARF file as copy of original ARF
    char arf_filepathname[MAXFILENAME];
		strcpy(arf_filepathname, inst->filepath);
		strcat(arf_filepathname, inst->tel->arf_filename);
    
    arf_out = initARFFile(&par, arf_filepathname, &status);
    CHECK_STATUS_BREAK(status);

    // --- End of Initialization ---


    // --- ARF generation ---
    headas_chat(2, "Calculating ARF correction ...\n");
    sixte_arfgen(inst, memelf, ac, arf_out, &par, &status);

    // --- End of ARF generation ---

  } while(0); // END of ERROR HANDLING Loop

  // --- Clean up ---
  headas_chat(2, "\ncleaning up ...\n");

  // Clean up the random number generator
  sixt_destroy_rng();

  // Release memory
  sixte_arfgen_freepar(&par);
  freeARFFile(&arf_out, &status);
  freeEventFile(&elf, &status);
  freeEventFile(&memelf, &status);
  freeAttitude(&ac);
  destroyGenInst(&inst, &status);

  // Remove temporary event file
  if (access(par.eventlist_filename, F_OK) == 0) {
    status = remove(par.eventlist_filename);
  }

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


////////////////////////////////////////////////////////////////////////////////
// Function definitions.
////////////////////////////////////////////////////////////////////////////////

void sixte_arfgen(GenInst* const inst, EventFile* elf,
            Attitude* const ac, ARFFile* arf_out,
            Parameters* par, int* status) {
  // Check for previous error
  CHECK_STATUS_VOID(*status);

  // Get pointer to ARF
  struct ARF* arf = inst->tel->arf;

  // Allocate memory for ARF correction factors
  ARFCorr* arf_corr = initARFCorr(arf, status);

  // Calculate ARF correction factors
  calc_arf_corr(arf_corr, arf, inst, ac, elf, par, status);

  // Update ARF effective areas with corrected values
  update_arf(arf_out, arf_corr, arf, status);

  // Free memory
  freeARFCorr(&arf_corr);
}


void calc_arf_corr(ARFCorr* arf_corr, const struct ARF* const arf,
                   GenInst* const inst, Attitude* const ac, EventFile* elf,
                   Parameters* par, int* status)
{
  // Allocate memory for row_status, required by fits filtering function
  char* row_status = malloc(par->n_photons * sizeof(*row_status));

  // Loop over specified ARF bins
  for (int ii = 0; ii < arf_corr->n_corr; ii++) {
    // Get energies of this ARF bin
    double e_low = arf->LowEnergy[arf_corr->bin_idxs[ii]];
    double e_hi = arf->HighEnergy[arf_corr->bin_idxs[ii]];

    // Simulate n_photons photons for this bin
    for (int ph_idx = 0; ph_idx < par->n_photons; ph_idx++) {
      // Generate photon within specified region and energy range
      Photon ph;
      double time = ph_idx*(1./par->photon_rate);
      gen_photon(&ph, e_low, e_hi,
                 par->SourceRA, par->SourceDec, time, status);

      // Photon imaging
      Impact imp;
      int isimg = phimg(inst->tel, ac, &ph, &imp, status);

      // If the photon is not imaged but lost in the optical system,
      // continue with the next one
      if (0 == isimg) continue;

      // Determine rawx/rawy [pixels]
      int rawx, rawy;
      double xr, yr;
      getGenDetAffectedPixel(inst->det->pixgrid, imp.position.x, imp.position.y,
                             &rawx, &rawy, &xr, &yr);

      // Add to temporary event file if valid pixel values
  		if ((rawx >= 0) && (rawy >= 0)) {
  			Event event = {.rawx = rawx, .rawy = rawy, .ph_id = {imp.ph_id},
                       .src_id = {imp.src_id}, .time = imp.time};
  			addEvent2File(elf, &event, status);
  		}

      CHECK_STATUS_VOID(*status);
    }

    // Run the event projection
    phproj(inst, ac, elf, par->TSTART, par->TSTART + par->exposure, status);

    // Add X,Y information to event file
    addXY2eventfile(elf, &(par->RefRA), &(par->RefDec), par->Projection, status);

    // Apply region filter
    long n_good_rows = -1;
    fits_find_rows(elf->fptr, par->regfilter, 1, elf->nrows, &n_good_rows,
                   row_status, status);
    assert( (n_good_rows >= 0) && (n_good_rows <= elf->nrows) );

    // Calculate correction factor for this energy bin
    arf_corr->corr_fac[ii] = (double)n_good_rows/par->n_photons;

    // Clear elf for next iteration
    fits_delete_rows(elf->fptr, 1, elf->nrows, status);
    elf->nrows = 0;

    // Print progress status
    headas_chat(2, "\r%.0f %%", 100.*(ii+1)/arf_corr->n_corr);
    fflush(NULL);
  }

  // Cleanup
  free(row_status);
}


void update_arf(ARFFile* arf_out, ARFCorr* arf_corr,
                const struct ARF* const arf, int* status) {
  // Allocate and initialize gsl interpolation object
  gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, arf_corr->n_corr);
  gsl_interp_init(interp, arf_corr->corr_energ, arf_corr->corr_fac,
                  arf_corr->n_corr);

  // Allocate and initialize gsl accelerator object
  gsl_interp_accel* acc = gsl_interp_accel_alloc ();

  // Load old specref (of the original ARF)
  double* specref = malloc(arf->NumberEnergyBins * sizeof(*specref));
  int anynul;
  fits_read_col(arf_out->fptr, TDOUBLE, arf_out->cspecresp, 1,
                1, arf->NumberEnergyBins, NULL, specref,
                &anynul, status);

  // Loop over all ARF bins
  for (int bin_idx = 0; bin_idx < arf->NumberEnergyBins; bin_idx++) {
    // Get energies of this ARF bin
    double e_low = arf->LowEnergy[bin_idx];
    double e_hi = arf->HighEnergy[bin_idx];

    // Calculate correction factor for this bin by interpolation
    double arf_corr_interp = gsl_interp_eval(interp, arf_corr->corr_energ,
                               arf_corr->corr_fac, 0.5*(e_low+e_hi), acc);

    // Calculate corrected effective area
    specref[bin_idx] *= arf_corr_interp;
  }

  // Update specref (i.e., effective area) values in output ARF file
  fits_write_col(arf_out->fptr, TDOUBLE, arf_out->cspecresp,
                 1, 1, arf->NumberEnergyBins, specref, status);

  // Cleanup
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);
  free(specref);
}


ARFFile* initARFFile(const Parameters* const par, char* arf_filepathname, int* status) {
  CHECK_STATUS_RET(*status, NULL);

  // Allocate memory
  ARFFile* arf_file = malloc(sizeof *arf_file);

  // Create the fits file
  fits_create_file_clobber(&(arf_file->fptr), par->ARFCorr, par->clobber, status);

  // Init column numbers
  arf_file->cenerg_lo = 1;
  arf_file->cenerg_hi = 2;
  arf_file->cspecresp = 3;

  // Copy data from original ARF. Specref will be updated later.
  fitsfile* fptr_orig;
  fits_open_table(&fptr_orig, arf_filepathname, READONLY, status);
  fits_copy_file(fptr_orig, arf_file->fptr, 1, 1, 1, status);
  fits_close_file(fptr_orig, status);
  fits_flush_file(arf_file->fptr, status);

  // Check if fits operations have been successful
  if (*status != EXIT_SUCCESS) {
    fits_report_error(stderr, *status);
    return NULL;
  }

  return arf_file;
}


void freeARFFile(ARFFile** arf_file, int* status) {
  if (*arf_file != NULL) {
    if ((*arf_file)->fptr != NULL) {
      fits_close_file_chksum((*arf_file)->fptr, status);
    }
    free(*arf_file);
    *arf_file = NULL;
  }
}


ARFCorr* initARFCorr(const struct ARF* const arf, int* status) {
  ARFCorr* arf_corr = malloc(sizeof(ARFCorr));
  CHECK_MALLOC_RET_NULL_STATUS(arf_corr, *status);

  // Since we only simulate every n-th ARF bin (where n is given by sampling_factor),
  // we need space for NumberEnergyBins/sampling_factor (rounded up) correction
  // factors to cover the whole ARF (rounding up to always including last bin and
  // avoid extrapolation when writing the corrected ARF).
  arf_corr->n_corr = ceil(1.*arf->NumberEnergyBins/sampling_factor);

  // Allocate memory for members.
  arf_corr->corr_fac = malloc(arf_corr->n_corr * sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(arf_corr->corr_fac, *status);

  arf_corr->corr_energ = malloc(arf_corr->n_corr * sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(arf_corr->corr_energ, *status);

  arf_corr->bin_idxs = malloc(arf_corr->n_corr * sizeof(int));
  CHECK_MALLOC_RET_NULL_STATUS(arf_corr->bin_idxs, *status);

  // Initialize bin_idxs and corr_energ.
  for (int ii = 0; ii < arf_corr->n_corr; ii++) {
    // Calculate bin index (at the moment, this is just every n-th bin)
    int bin_idx;
    if (ii < arf_corr->n_corr-1) {
      // Simulate every n-th bin (where n is given by sampling_factor)
      bin_idx = ii*sampling_factor;
    } else {
      // Use last ARF bin in the last iteration
      bin_idx = arf->NumberEnergyBins-1;
    }
    arf_corr->bin_idxs[ii] = bin_idx;

    // Save energy in corr_energ (midpoint between e_low and e_hi)
    arf_corr->corr_energ[ii] = 0.5 * (arf->LowEnergy[bin_idx]
                                      + arf->HighEnergy[bin_idx]);
  }

  return arf_corr;
}


void freeARFCorr(ARFCorr** arf_corr) {
  if (*arf_corr != NULL) {
    free((*arf_corr)->bin_idxs);
    free((*arf_corr)->corr_fac);
    free((*arf_corr)->corr_energ);
  }
  free(*arf_corr);
  *arf_corr = NULL;
}


void gen_photon(Photon* ph, float E_low, float E_high, double RA,
                double Dec, double time, int* status) {
  // Set photon position
  set_ph_position(ph, RA, Dec);

  // Set photon energy
  set_ph_energy(ph, E_low, E_high, status);

  // Set additional photon information (time, ph_id, src_id)
  set_ph_info(ph, time);
}


void set_ph_info(Photon* ph, double time) {
  // Set photon time
  ph->time = time;

  // Set photon id
  static long ph_id = 0;
  ph->ph_id = ph_id;
  ph_id++;

  // Set src_id
  ph->src_id = 0;
}


void set_ph_energy(Photon* ph, double E_low, double E_high, int* status) {
  // Sample an energy uniformly in [E_low, E_high)
  double energy = (E_high - E_low)*sixt_get_random_number(status) + E_low;

  // Set the photon energy
  ph->energy = energy;
}


void set_ph_position(Photon* ph, double RA, double Dec) {
  ph->ra = RA*M_PI/180.;
  ph->dec = Dec*M_PI/180.;
  return;
}


// Reads WCS information from an image
static void getWCSfromImage(char* ImageFile, char* Projection, float* RefRA,
                            float* RefDec, int* status) {
  // Open ImageFile
  fitsfile* fptr;
  fits_open_file(&fptr, ImageFile, READONLY, status);

  // Read WCS keywords
  *RefRA = 5;
  *RefDec = 6;
  fits_read_key(fptr, TFLOAT, "CRVAL1", RefRA, NULL, status);
  fits_read_key(fptr, TFLOAT, "CRVAL2", RefDec, NULL, status);

  char ctype1[MAXMSG];
  fits_read_key(fptr, TSTRING, "CTYPE1", ctype1, NULL, status);
  if (strncmp(ctype1, "RA-", 3) != 0) {
    SIXT_ERROR("CTYPE1 does not start with \"RA\"");
  }

  // Find projection type (e.g., as in RA---TAN or RA---AIT)
  int idx = 2; // Start at first -
  while(ctype1[idx] == 45) idx++; // (45 is ASCII encoding for -)
  strcpy(Projection, &ctype1[idx]);
  assert(strlen(Projection) == 3);
}

int sixte_arfgen_getpar(Parameters* const par) {
  // Error status
  int status = EXIT_SUCCESS;

  // Get parameters
  query_simput_parameter_int("Seed", &(par->Seed), &status);
  query_simput_parameter_bool("clobber", &(par->clobber), &status);
  query_simput_parameter_file_name("ARFCorr", &(par->ARFCorr), &status);
  query_simput_parameter_file_name("XMLFile", &(par->XMLFile), &status);
  query_simput_parameter_double("MJDREF", &(par->MJDREF), &status);
  query_simput_parameter_int("n_photons", &(par->n_photons), &status);
  query_simput_parameter_double("photon_rate", &(par->photon_rate), &status);

  query_simput_parameter_file_name("Attitude", &(par->Attitude), &status);

  // only load PointingRA, PointingDec if Attitude is not given
  if (par->Attitude) {
    // set to default values
    par->PointingRA=0.0;
    par->PointingDec=0.0;
    par->rollangle=0.0;
    headas_chat(2, "using Attitude File: %s \n", par->Attitude);

    query_simput_parameter_double("TSTART", &(par->TSTART), &status);
    query_simput_parameter_double("Exposure", &(par->exposure), &status);

    // If exposure is set, make sure TSTART is also set (and vice versa)
    if ( ((par->exposure >= 0) && (par->TSTART < 0)) ||
         ((par->exposure < 0) && (par->TSTART >= 0)) ) {
      SIXT_ERROR("Must provide TSTART AND Exposure");
      status = EXIT_FAILURE;
      return status;
    }
  } else {
    query_simput_parameter_double("PointingRA", &(par->PointingRA), &status);
    query_simput_parameter_double("PointingDec", &(par->PointingDec), &status);
    query_simput_parameter_double("rollangle", &(par->rollangle), &status);
  }

  query_simput_parameter_double("SourceRA", &(par->SourceRA), &status);
  query_simput_parameter_double("SourceDec", &(par->SourceDec), &status);

  query_simput_parameter_string_buffer("Projection", par->Projection, MAXFILENAME,
                                       &status);
  if ( strncasecmp(par->Projection, "none", 5) == 0 ) {
    // No projection specified, so read WCS keywords from image
    query_simput_parameter_file_name("ImageFile", &(par->ImageFile), &status);
    getWCSfromImage(par->ImageFile, par->Projection, &par->RefRA, &par->RefDec,
                    &status);
  } else {
    query_simput_parameter_float("RefRa", &par->RefRA, &status);
    query_simput_parameter_float("RefDec", &par->RefDec, &status);
  }

  char* buf;
  query_simput_parameter_file_name("regfile", &buf, &status);
  // Save as regfilter, using appropriate syntax for fits row filtering
  sprintf(par->regfilter, "regfilter('%s')", buf);
  free(buf);

  strcpy(par->eventlist_filename, "tmp_evt.fits");

  return(status);
}


void sixte_arfgen_freepar(Parameters* par) {
  free(par->ARFCorr);
  free(par->XMLFile);
  free(par->ImageFile);
}
