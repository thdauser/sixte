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

#include "phoimg.h"


////////////////////////////////////
/** Main procedure. */
int phoimg_main() {
  struct Parameters par;

  Attitude* ac=NULL;

  PhotonFile* plf=NULL;
  ImpactFile* ilf=NULL;

  // Instrument data structure including telescope information like the PSF,
  // vignetting function, focal length, and FOV diameter.
  GenInst* inst=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("phoimg");
  set_toolversion("0.04");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=phoimg_getpar(&par))) break;

    headas_chat(3, "initialize ...\n");

    // Determine the photon list file name.
    char photonlist_filename[MAXFILENAME];
    strcpy(photonlist_filename, par.PhotonList);

    // Determine the impact list output file.
    char impactlist_filename[MAXFILENAME];
    strcpy(impactlist_filename, par.ImpactList);

    // Initialize the random number generator.
    unsigned int seed=getSeed(par.Seed);
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the instrument XML definition file.
    char xml_filename[MAXFILENAME];
    strcpy(xml_filename, par.XMLFile);

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
			     par.RA*M_PI/180., par.Dec*M_PI/180., par.rollangle*M_PI/180.,&status);
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

    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(3, "start imaging process ...\n");

    // Open the input photon list file.
    plf=openPhotonFile(photonlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read header keywords.
    char telescop[MAXMSG], instrume[MAXMSG], comment[MAXMSG];
    fits_read_key(plf->fptr, TSTRING, "TELESCOP", &telescop, comment, &status);
    fits_read_key(plf->fptr, TSTRING, "INSTRUME", &instrume, comment, &status);
    CHECK_STATUS_BREAK(status);

    double mjdref, timezero, tstart, tstop;
    fits_read_key(plf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(plf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(plf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(plf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output impact list file.
    ilf=openNewImpactFile(impactlist_filename,
			  telescop, instrume,
			  inst->tel->arf->Filter,
			  inst->tel->arf_filename,
			  inst->det->rmf_filename,
			  mjdref, timezero, tstart, tstop,
			  par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);
    CHECK_STATUS_BREAK(status);

    // Scan the entire photon list.
    int progress=0;
    while (plf->row < plf->nrows) {

      Photon photon={.time=0.};

      // Read an entry from the photon list:
      status=PhotonFile_getNextRow(plf, &photon);
      CHECK_STATUS_BREAK(status);

      // Check whether we are still within the requested time interval.
      if (photon.time < par.TSTART) continue;
      if (photon.time > par.TSTART+par.Exposure) break;

      // Photon Imaging.
      Impact impact;
      int isimg=phimg(inst->tel, ac, &photon, &impact, &status);
      CHECK_STATUS_BREAK(status);

      // If the photon is not imaged but lost in the optical system,
      // continue with the next one.
      if (0==isimg) continue;

      // Write the impact to the output file.
      addImpact2File(ilf, &impact, &status);
      CHECK_STATUS_BREAK(status);

      // Program progress output.
      while ((int)((photon.time-par.TSTART)*1000./par.Exposure)>progress) {
	progress++;
	headas_chat(2, "\r%.1lf %%", progress*1./10.);
	fflush(NULL);
      }

    }; // END of photon processing loop.
    CHECK_STATUS_BREAK(status);

    // Progress output.
    headas_chat(2, "\r%.1lf %%\n", 100.);
    fflush(NULL);

    // --- END of imaging process ---

  } while(0); // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Clean up the random number generator.
  sixt_destroy_rng();

  // Close the FITS files.
  freeImpactFile(&ilf, &status);
  freePhotonFile(&plf, &status);

  freeAttitude(&ac);
  destroyGenInst(&inst, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int phoimg_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the photon list");
    return(status);
  }
  strcpy(par->PhotonList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("ImpactList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the impact list");
    return(status);
  }
  strcpy(par->ImpactList, sbuffer);
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
    SIXT_ERROR("failed reading the name of the attitude file");
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

  query_simput_parameter_float("rollangle",&(par->rollangle),&status);


  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading MJDREF");
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
