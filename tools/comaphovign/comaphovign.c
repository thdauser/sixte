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

#include "comaphovign.h"


////////////////////////////////////
/** Main procedure. */
int comaphovign_main() {
  struct Parameters par;

  Attitude* ac=NULL;
  PhotonFile* plif=NULL;

  fitsfile* ofptr=NULL;
  long onrows=0;
  int cenergy, cra, cdec, ctime;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("comaphovign");
  set_toolversion("0.07");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=comaphovign_getpar(&par))) break;

    headas_chat(3, "initialize ...\n");

    // Determine the input photon list file name.
    char inputlist_filename[MAXFILENAME];
    strcpy(inputlist_filename, par.InputList);

    // Determine the output photon list file name.
    char outputlist_filename[MAXFILENAME];
    strcpy(outputlist_filename, par.OutputList);

    // Determine the random number seed.
    int seed=getSeed(par.Seed);

    // Initialize the random number generator.
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Open the input photon list file.
    plif=openPhotonFile(inputlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read keywords.
    char comment[MAXMSG];
    double mjdref=0.;
    fits_read_key(plif->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'MJDREF' from input "
	      "photon list '%s'", inputlist_filename);
      SIXT_ERROR(msg);
      break;
    }

    double timezero=0.;
    fits_write_errmark();
    int status2=EXIT_SUCCESS;
    fits_read_key(plif->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status2);
    fits_clear_errmark();
    if (EXIT_SUCCESS!=status2) {
      timezero=0.;
    }

    double tstart=0.;
    fits_read_key(plif->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTART' from input "
	      "photon list '%s'", inputlist_filename);
      SIXT_ERROR(msg);
      break;
    }

    double tstop=0.;
    fits_read_key(plif->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTOP' from input "
	      "photon list '%s'", inputlist_filename);
      SIXT_ERROR(msg);
      break;
    }

    // Make sure that TIMEZERO==0.0.
    verifyTIMEZERO(timezero, &status);
    CHECK_STATUS_BREAK(status);

    // Set up the Attitude.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a simple pointing attitude.
      ac=getPointingAttitude(mjdref, tstart, tstop,
			     par.RA*M_PI/180., par.Dec*M_PI/180., &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval is a subset of the period
      // described by the attitude file.
      checkAttitudeTimeCoverage(ac, mjdref, tstart, tstop, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of setting up the attitude.

    // Open the output photon list file.
    // Check if the file already exists.
    int exists;
    fits_file_exists(outputlist_filename, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(outputlist_filename);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", outputlist_filename);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }

    // Create a new empty photon list FITS file.
    fits_create_file(&ofptr, outputlist_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Create the binary table for the photon list.
    char *ttype[] = { "ENERGY", "RA", "DEC" };
    char *tform[] = { "E", "E", "E" };
    char *tunit[] = { "keV", "deg", "deg" };
    fits_create_tbl(ofptr, BINARY_TBL, 0, 3, ttype, tform, tunit, 
		    "PHOTONS", &status);
    CHECK_STATUS_BREAK(status);

    // Check if the time information should be stored in the 
    // output file.
    if (0!=par.TimeColumn) {
      fits_insert_col(ofptr, 1, "TIME", "D", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Determine the column numbers in the output file.
    fits_get_colnum(ofptr, CASEINSEN, "ENERGY", &cenergy, &status); 
    fits_get_colnum(ofptr, CASEINSEN, "RA", &cra, &status);
    fits_get_colnum(ofptr, CASEINSEN, "DEC", &cdec, &status);
    if (0!=par.TimeColumn) {
      fits_get_colnum(ofptr, CASEINSEN, "TIME", &ctime, &status);
    }
    CHECK_STATUS_BREAK(status);

    // Add header information about program parameters.
    // The second parameter "1" means that the headers are written
    // to the first extension.
    HDpar_stamp(ofptr, 1, &status);
    CHECK_STATUS_BREAK(status);

    // Move back to the second HDU.
    fits_movabs_hdu(ofptr, 2,NULL, &status);
    CHECK_STATUS_BREAK(status);
    
    // Write FITS header keywords.
    fits_update_key(ofptr, TSTRING, "ATTITUDE",
		    par.Attitude, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(ofptr, TDOUBLE, "MJDREF",
		    &mjdref, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(ofptr, TDOUBLE, "TIMEZERO", 
		    &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---

    headas_chat(3, "apply vignetting ...\n");

    // Scan the entire photon list.
    int progress=0;  
    while (plif->row < plif->nrows) {

      Photon photon={.time=0.};
      
      // Read an entry from the photon list:
      status=PhotonFile_getNextRow(plif, &photon);
      CHECK_STATUS_BREAK(status);

      // Apply the vignetting.
      // Compare the photon direction to the direction of the telescope axis.
      Vector nz=getTelescopeNz(ac, photon.time, &status);
      CHECK_STATUS_BREAK(status);
      Vector phodir=unit_vector(photon.ra, photon.dec);

      double cosine=nz.x*phodir.x + nz.y*phodir.y + nz.z*phodir.z;
      double rnd=sixt_get_random_number(&status);
      CHECK_STATUS_BREAK(status);

      // Delete the photon from the FITS file.
      if (cosine>=rnd) {
	// Append the photon to the output file.

	// Insert a new, empty row to the table:
	fits_insert_rows(ofptr, onrows, 1, &status);
	CHECK_STATUS_BREAK(status);
	onrows++;
	
	// Store the data in the FITS file.
	if (0!=par.TimeColumn) {
	  fits_write_col(ofptr, TDOUBLE, ctime, 
			 onrows, 1, 1, &photon.time, &status);
	}
	fits_write_col(ofptr, TFLOAT, cenergy, 
		       onrows, 1, 1, &photon.energy, &status);
	float fbuffer=photon.ra * 180./M_PI;
	fits_write_col(ofptr, TFLOAT, cra, 
		       onrows, 1, 1, &fbuffer, &status);
	fbuffer=photon.dec * 180./M_PI;
	fits_write_col(ofptr, TFLOAT, cdec, 
		       onrows, 1, 1, &fbuffer, &status);
	CHECK_STATUS_BREAK(status);
      }

      // Program progress output.
      while ((int)(plif->row*1000./plif->nrows)>progress) {
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
  if (NULL!=ofptr) {
    fits_close_file(ofptr, &status);
    ofptr=NULL;
  }
  freePhotonFile(&plif, &status);
  freeAttitude(&ac);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int comaphovign_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("InputList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("could not read the name of the input photon list");
    return(status);
  } 
  strcpy(par->InputList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("OutputList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("could not read the name of the output photon list");
    return(status);
  } 
  strcpy(par->OutputList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("could not read the name of the attitude");
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

  status=ape_trad_query_bool("TimeColumn", &par->TimeColumn);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the TimeColumn parameter");
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

