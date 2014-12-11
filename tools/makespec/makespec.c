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

#include "makespec.h"


////////////////////////////////////
/** Main procedure. */
int makespec_main() {
  // Program parameters.
  struct Parameters par; 
  
  // Full file name with filter.
  char evtlistfiltered[2*MAXFILENAME];

  // Input event file.
  fitsfile* ef=NULL;

  // Output spectrum.
  long* spec=NULL;
  fitsfile* sf=NULL;

  // GTI.
  GTI* gti=NULL;

  // Instrument response.
  struct RMF* rmf=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("makespec");
  set_toolversion("0.09");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=makespec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");
    
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.EventFilter);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(evtlistfiltered, par.EventList);
    }else{
      if(0>=sprintf(evtlistfiltered, "%s[EVENTS][%s]", par.EventList, par.EventFilter)){
	status=EXIT_FAILURE;
	SIXT_ERROR("Assembling file name failed.");
	break;
      }
    }

    // Open the input event file.
    fits_open_table(&ef, evtlistfiltered, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read required keywords.
    char comment[MAXMSG];
    char telescop[MAXMSG];
    fits_read_key(ef, TSTRING, "TELESCOP", telescop, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'TELESCOP' in event file");
      break;
    }
    char instrume[MAXMSG];
    fits_read_key(ef, TSTRING, "INSTRUME", instrume, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'INSTRUME' in event file");
      break;
    }
    char filter[MAXMSG];
    fits_read_key(ef, TSTRING, "FILTER", filter, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'FILTER' in event file");
      break;
    }

    char ancrfile[MAXMSG];
    fits_read_key(ef, TSTRING, "ANCRFILE", ancrfile, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'ANCRFILE' in event file");
      break;
    }
    char respfile[MAXMSG];
    fits_read_key(ef, TSTRING, "RESPFILE", respfile, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'RESPFILE' in event file");
      break;
    }

    char date_obs[MAXMSG];
    fits_read_key(ef, TSTRING, "DATE-OBS", date_obs, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'DATE-OBS' in event file");
      break;
    }
    char time_obs[MAXMSG];
    fits_read_key(ef, TSTRING, "TIME-OBS", time_obs, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'TIME-OBS' in event file");
      break;
    }
    char date_end[MAXMSG];
    fits_read_key(ef, TSTRING, "DATE-END", date_end, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'DATE-END' in event file");
      break;
    }
    char time_end[MAXMSG];
    fits_read_key(ef, TSTRING, "TIME-END", time_end, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not find keyword 'TIME-END' in event file");
      break;
    }

    // Load the GTI extension in order to be able to determine the 
    // exposure time.
    gti=loadGTI(par.EventList, &status);
    CHECK_STATUS_BREAK(status);
    double exposure=sumGTI(gti);

    // Determine the column containing the signal information.
    int csignal;
    fits_get_colnum(ef, CASEINSEN, "SIGNAL", &csignal, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the number of rows.
    long nrows;
    fits_get_num_rows(ef, &nrows, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the random number generator seed.
    int seed;
    if (-1!=par.Seed) {
      seed=par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed=(int)time(NULL);
    }

    // Load the EBOUNDS of the RMF.
    char filepathname[MAXFILENAME];
    if (strlen(par.RSPPath)>0) {
      strcpy(filepathname, par.RSPPath);
      strcat(filepathname, "/");
      strcat(filepathname, respfile);
    } else {
      // The file should be located in the working directory.
      strcpy(filepathname, respfile);
    }
    struct RMF* rmf=getRMF(&status);
    CHECK_STATUS_BREAK(status);
    loadEbounds(rmf, filepathname, &status);
    CHECK_STATUS_BREAK(status);

    // Initialize the random number generator.
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Allocate memory for the output spectrum.
    headas_chat(5, "create empty spectrum with %ld channels ...\n",
		rmf->NumberChannels);
    spec=(long*)malloc(rmf->NumberChannels*sizeof(long));
    CHECK_NULL_BREAK(spec, status, "memory allocation for spectrum failed");

    // Initialize the spectrum with 0.
    long ii; 
    for (ii=0; ii<rmf->NumberChannels; ii++) {
      spec[ii]=0;
    }

    // --- END of Initialization ---


    // --- Begin Spectrum Binning ---
    headas_chat(3, "calculate spectrum ...\n");

    // LOOP over all events in the FITS table.
    for (ii=0; ii<nrows; ii++) {
      
      // Read the next event from the file.
      float signal;
      float fnull=0.0;
      int anynul=0;
      fits_read_col(ef, TFLOAT, csignal, ii+1, 1, 1, 
		    &fnull, &signal, &anynul, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the PHA channel.
      long pha=getEBOUNDSChannel(signal, rmf);
      
      // Add the event to the spectrum.
      long idx=pha-rmf->FirstChannel;
      if(idx>=0) {
	assert(idx<rmf->NumberChannels);      
	spec[idx]++;
      }
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

    // Store the spectrum in the output file.
    headas_chat(3, "store spectrum ...\n");

    // Create a new FITS-file (remove existing one before):
    remove(par.Spectrum);
    char buffer[MAXFILENAME];
    sprintf(buffer, "%s(%s%s)", par.Spectrum, SIXT_DATA_PATH, 
	    "/templates/makespec.tpl");
    fits_create_file(&sf, buffer, &status);
    CHECK_STATUS_BREAK(status);

    // Move to the HDU containing the binary table.
    int hdutype;
    fits_movabs_hdu(sf, 2, &hdutype, &status);
    CHECK_STATUS_BREAK(status);

    // Get column numbers.
    int cchannel, ccounts;
    fits_get_colnum(sf, CASEINSEN, "CHANNEL", &cchannel, &status);
    fits_get_colnum(sf, CASEINSEN, "COUNTS", &ccounts, &status);
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_update_key(sf, TSTRING, "TELESCOP", telescop, "", &status);
    fits_update_key(sf, TSTRING, "INSTRUME", instrume, "", &status);
    fits_update_key(sf, TSTRING, "FILTER", filter, "", &status);
    fits_update_key(sf, TSTRING, "DATE-OBS", date_obs, "", &status);
    fits_update_key(sf, TSTRING, "TIME-OBS", time_obs, "", &status);
    fits_update_key(sf, TSTRING, "DATE-END", date_end, "", &status);
    fits_update_key(sf, TSTRING, "TIME-END", time_end, "", &status);
    fits_update_key(sf, TSTRING, "ANCRFILE", ancrfile, 
		    "ancillary response file", &status);
    fits_update_key(sf, TSTRING, "RESPFILE", respfile, 
		    "response file", &status);
    fits_update_key(sf, TSTRING, "BACKFILE", "", 
		    "background file", &status);
    fits_update_key(sf, TLONG, "DETCHANS", &rmf->NumberChannels,
		    "number of detector channels", &status);
    fits_update_key(sf, TSTRING, "CORRFILE", "", 
		    "none", &status);
    fits_update_key(sf, TDOUBLE, "EXPOSURE", &exposure,
		    "exposure time", &status);
    fits_update_key(sf, TSTRING, "FilterExpr", par.EventFilter,
		    "used filter expression", &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all channels in the spectrum.
    for (ii=0; ii<rmf->NumberChannels; ii++) {    
      long channel=ii+rmf->FirstChannel;
      fits_write_col(sf, TLONG, cchannel, ii+1, 1, 1, &channel, &status);
    }
    fits_write_col(sf, TLONG, ccounts, 1, 1, rmf->NumberChannels,
		   spec, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  if (NULL!=ef) fits_close_file(ef, &status);
  if (NULL!=sf) fits_close_file(sf, &status);

  // Release memory.
  if (NULL!=spec) free(spec);
  freeRMF(rmf);
  freeGTI(&gti);

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int makespec_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list file");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("EventFilter", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the event filter expression");
    return(status);
  } 
  strcpy(par->EventFilter, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Spectrum", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output spectrum file");
    return(status);
  } 
  strcpy(par->Spectrum, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("RSPPath", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading path to the response files");
    return(status);
  } 
  strcpy(par->RSPPath, sbuffer);
  free(sbuffer);

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
