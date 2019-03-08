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


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "htrssim.h"

////////////////////////////////////
/** Main procedure. */
int htrssim_main() {

  // Containing all programm parameters read by PIL
  struct Parameters parameters;

  // HTRSDetector data structure.
  // This is the main object for the simulation of the HTRS. It contains all
  // important parameters of the detector and is used to create an event file
  // from the incident photon events.
  // Before the first usage it has to be initialized.
  HTRSDetector detector;

  ImpactFile* impactlistfile=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("htrssim");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;

    // Initialize the random number generator.
    sixt_init_rng((int)time(NULL), &status);
    CHECK_STATUS_BREAK(status);

    // Open the impact list FITS file.
    impactlistfile=openImpactFile(parameters.impactlist_filename,
				  READONLY, &status);
    CHECK_STATUS_BREAK(status);


    // Detector settings.
    // Store the settings for the HTRSDetector in the corresponding
    // data structure and call the initialization routine in order
    // to allocate memory etc.
#ifdef HTRS_HEXPIXELS
    struct HTRSDetectorParameters hdparameters = {
      .pixels = { .npixels = 37,
		  .pixelwidth = parameters.pixelwidth },
      .generic = { .ccsigma = parameters.ccsigma,
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .arf_filename = parameters.arf_filename,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .slow_shaping_time  = parameters.slow_shaping_time,
      .fast_shaping_time  = parameters.fast_shaping_time,
      .reset_time         = parameters.reset_time,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template /* String address!! */
    };
#endif
#ifdef HTRS_ARCPIXELS
    int nrings = 4;
    // Out-of-focus distance 12 cm (detector radius 14.15 mm)

    // Configuration with 31 pixels with each pixel having the same area.
    int npixels[4] = { 1, 6, 12, 12 };
    double radii[4] = { 2.16e-3, 5.70e-3, 9.39e-3, 12.0e-3 };
    double offset_angles[4] = { 0., 0., 0., 0. };

    struct HTRSDetectorParameters hdparameters = {
      .pixels = { .nrings = nrings,
		  .npixels = npixels,
		  .radius = radii,
		  .offset_angle = offset_angles,
		  .mask_spoke_width = parameters.mask_spoke_width },
      .generic = { .ccsigma = parameters.ccsigma,
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .arf_filename = parameters.arf_filename,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .slow_shaping_time  = parameters.slow_shaping_time,
      .fast_shaping_time  = parameters.fast_shaping_time,
      .reset_time         = parameters.reset_time,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template /* String address!! */
    };
#endif
    if(EXIT_SUCCESS!=(status=initHTRSDetector(&detector, &hdparameters))) break;
    // END of DETECTOR CONFIGURATION SETUP

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");
    Impact impact; // Buffer to store the impacts read from the FITS file.

    // Loop over all impacts in the FITS file.
    while (impactlistfile->row<impactlistfile->nrows) {

      getNextImpactFromFile(impactlistfile, &impact, &status);
      if (EXIT_SUCCESS!=status) break;

      // Check whether the event lies in the specified time interval:
      if ((impact.time<parameters.TSTART)||
	  (impact.time>parameters.TSTART+parameters.timespan)) continue;

      // Call the photon detection routine that generates the right charge
      // and stores it in the detector pixels.
      // Before generating and adding the charge to the detector the routine also
      // checks, whether the integration time is exceeded and performs the readout
      // in that case.
      status=addImpact2HTRSDetector(&detector, &impact);
      if (EXIT_SUCCESS!=status) break;

    } // END of scanning the impact list.
    if (EXIT_SUCCESS!=status) break;


    // Assign event grades to the detected events.
    status=HTRSassignEventGrades(detector);
    if (EXIT_SUCCESS!=status) break;

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Close the FITS files.
  freeImpactFile(&impactlistfile, &status);

  // Release memory of detector.
  status+=cleanupHTRSDetector(&detector);

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
static int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the name of the impact list file (FITS file).
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
  }

  // Get the slow shaping time for the pulses.
  else if ((status = PILGetReal("slow_shaping_time", &parameters->slow_shaping_time))) {
    HD_ERROR_THROW("Error reading the slow shaping time for pulses!\n", status);
  }

  // Get the fast shaping time for the pulses.
  else if ((status = PILGetReal("fast_shaping_time", &parameters->fast_shaping_time))) {
    HD_ERROR_THROW("Error reading the fast shaping time for pulses!\n", status);
  }

  // Get the reset time for a detector pixel.
  else if ((status = PILGetReal("reset_time", &parameters->reset_time))) {
    HD_ERROR_THROW("Error reading the reset time for the detector pixels!\n", status);
  }

#ifdef HTRS_HEXPIXELS
  // [m]
  else if ((status = PILGetReal("pixelwidth", &parameters->pixelwidth))) {
    HD_ERROR_THROW("Error reading the width of the detector pixels!\n", status);
  }
#endif

  // [m]
  else if ((status = PILGetReal("ccsigma", &parameters->ccsigma))) {
    HD_ERROR_THROW("Error reading the charge cloud sigma!\n", status);
  }
  if (status) return(status);


  // Read the detector thresholds (either integer PHA or float energy):
  int pha_threshold;
  if ((status = PILGetInt("pha_threshold", &pha_threshold))) {
    HD_ERROR_THROW("Error: could not determine detector PHA threshold!\n", status);
    return(status);
  } else {
    parameters->pha_threshold = (long)pha_threshold;
  }
  if (parameters->pha_threshold==-1) {
    if ((status = PILGetReal4("energy_threshold", &parameters->energy_threshold))) {
      HD_ERROR_THROW("Error: could not determine detector energy threshold!\n", status);
      return(status);
    }
  } else {
    parameters->energy_threshold=0.;
  }

  // Get the name of the detector redistribution file (FITS file)
  if ((status = PILGetFname("rmf_filename", parameters->rmf_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector"
		   "redistribution matrix file (RMF)!\n", status);
  }

  if ((status = PILGetFname("arf_filename", parameters->arf_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector ARF!\n", status);
  }

  // Get the background count rate
  else if ((status = PILGetReal4("background_rate", &parameters->background_rate))) {
    HD_ERROR_THROW("Error: could not determine the detector background rate!\n", status);
  }

  // Get the name of the output event list (FITS file)
  else if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the event list file!\n", status);
  }

  // Get the start time of the simulation
  else if ((status = PILGetReal("TSTART", &parameters->TSTART))) {
    HD_ERROR_THROW("Error reading the 'TSTART' parameter!\n", status);
  }

  // Get the timespan for the simulation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

#ifdef HTRS_ARCPIXELS
  // Determine the width of the spokes of the HTRS mask.
  else if ((status = PILGetReal("mask_spoke_width", &parameters->mask_spoke_width))) {
    HD_ERROR_THROW("Error reading the spoke width of the HTRS mask!\n", status);
  }
#endif

  if (EXIT_SUCCESS!=status) return(status);

  // Set the event list template file:
  strcpy(parameters->eventlist_template, SIXT_DATA_PATH);
  strcat(parameters->eventlist_template, "/templates/htrs.eventlist.tpl");

  return(status);
}
