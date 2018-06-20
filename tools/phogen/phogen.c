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
*/

#include "phogen.h"


int phogen_main() 
{
  const double timezero=0.0;

  // Program parameters.
  struct Parameters par;

  // Instrument setup.
  GenInst* inst=NULL;

  // Attitude.
  Attitude* ac=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat=NULL;

  // Photon list file.
  PhotonFile* plf=NULL;
  
  // Error status.
  int status=EXIT_SUCCESS;  

  // Register HEATOOL.
  set_toolname("phogen");
  set_toolversion("0.05");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=phogen_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Determine the photon list output file.
    char photonlist_filename[MAXFILENAME];
    strcpy(photonlist_filename, par.PhotonList);

    // Initialize the random number generator.
    unsigned int seed=getSeed(par.Seed);
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile, 
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, seed, &status);
    CHECK_STATUS_BREAK(status);
    
    // Set up the Attitude.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a simple pointing attitude.
      ac=getPointingAttitude(par.MJDREF, par.TSTART, par.TSTART+par.Exposure,
			     par.RA*M_PI/180., par.Dec*M_PI/180., &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval for the simulation
      // is a subset of the period covered by the attitude file.
      checkAttitudeTimeCoverage(ac, par.MJDREF, par.TSTART, 
				par.TSTART+par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of setting up the attitude.

    // Load the SIMPUT X-ray source catalog.
    srccat=loadSourceCatalog(par.Simput, inst->tel->arf, &status);
    CHECK_STATUS_BREAK(status);

    // --- End of Initialization ---


    // --- Photon Generation Process ---

    // Open the output photon list file.
    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    plf=openNewPhotonFile(photonlist_filename,
			  telescop, instrume,
			  inst->tel->arf->Filter,
			  inst->tel->arf_filename,
			  inst->det->rmf_filename,
			  par.MJDREF, timezero, par.TSTART, par.TSTART+par.Exposure,
			  par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);

    // Start the actual photon generation (after loading required data):
    headas_chat(3, "start photon generation process ...\n");

    // Loop over photon generation, till the time of the photon exceeds 
    // the requested exposure time.
    // Simulation progress status (running from 0 to 1000).
    int progress=0;
    do {

      // Photon generation.
      Photon ph;
      int isph=phgen(ac, &srccat, 1, 
		     par.TSTART, par.TSTART+par.Exposure, 
		     par.MJDREF, par.dt, 
		     inst->tel->fov_diameter, &ph, &status);
      CHECK_STATUS_BREAK(status);

      // If no photon has been generated, break the loop.
      if (0==isph) break;

      // Check if the photon still is within the requested exposre time.
      if (ph.time>par.TSTART+par.Exposure) break;

      // Write the photon to the output file.
      status=addPhoton2File(plf, &ph);
      CHECK_STATUS_BREAK(status);

      // Program progress output.
      while ((int)((ph.time-par.TSTART)*1000./par.Exposure)>progress) {
	progress++;
	headas_chat(2, "\r%.1lf %%", progress*1./10.);
	fflush(NULL);
      }

    } while(1); // END of photon generation loop.
    CHECK_STATUS_BREAK(status);
 
    // Progress output.
    headas_chat(2, "\r%.1lf %%\n", 100.);
    fflush(NULL);

    // --- End of photon generation ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freePhotonFile(&plf, &status);
  freeSourceCatalog(&srccat, &status);
  freeAttitude(&ac);
  destroyGenInst(&inst, &status);

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}



int phogen_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the photon list");
    return(status);
  } 
  strcpy(par->PhotonList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mission", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the mission");
    return(status);
  } 
  strcpy(par->Mission, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Instrument", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the instrument");
    return(status);
  } 
  strcpy(par->Instrument, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mode", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the instrument mode");
    return(status);
  } 
  strcpy(par->Mode, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the attitude");
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the right ascension of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the declination of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the SIMPUT file");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading MJDREF");
    return(status);
  } 

  status=ape_trad_query_double("dt", &par->dt);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading dt");
    return(status);
  } 

  status=ape_trad_query_double("TSTART", &par->TSTART);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading TSTART");
    return(status);
  } 

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the exposure time");
    return(status);
  } 

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the seed for the random number generator");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}


