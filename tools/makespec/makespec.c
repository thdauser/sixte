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
  set_toolversion("0.16");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=makespec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Assemble event file name with filter
    if (0==strcasecmp(par.EventFilter,"NONE")) {
      strcpy(evtlistfiltered, par.EvtFile);
    }else{
      if(0>=sprintf(evtlistfiltered, "%s[EVENTS][%s]", par.EvtFile, par.EventFilter)){
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
    /* Try to obtain optional key "PIRMF" to be used for PI values */
    char pirmf[MAXMSG];
    fits_read_key(ef, TSTRING, "PIRMF", pirmf, comment, &status);
    if( status==VALUE_UNDEFINED || status==COL_NOT_FOUND || status==KEY_NO_EXIST ){
    	strcpy(pirmf,"NONE");
    	fits_clear_errmark();
    	status = EXIT_SUCCESS;
    }
    /* Try to obtain optional key "SPECARF" being the calibrated ARF */
    char specarf[MAXMSG];
    fits_read_key(ef, TSTRING, "SPECARF", specarf, comment, &status);
    if( status==VALUE_UNDEFINED || status==COL_NOT_FOUND || status==KEY_NO_EXIST ){
    	strcpy(specarf,"NONE");
    	fits_clear_errmark();
    	status = EXIT_SUCCESS;
    }

    // Load the GTI extension in order to be able to determine the
    // exposure time.
    gti=loadGTI(par.EvtFile, &status);
    CHECK_STATUS_BREAK(status);
    double exposure=sumGTI(gti);


    // Determine the column containing the signal information.
    char* pha2pi=NULL;

    int csignal;
    int coltmp;
    int usesignal = 0;

    char rmffile[MAXMSG];
    strcpy(rmffile, respfile );

    fits_get_colnum(ef, CASEINSEN, "PHA", &csignal, &status);
    // Print a warning if we find a PHA column in an X-IFU event file.
    if ( status == EXIT_SUCCESS && !strcmp(instrume, "XIFU") ) {
      SIXT_WARNING("Event file contains a PHA column and is not compatible with this version of makespec");
    }
    if( par.usepha==-1 || status == COL_NOT_FOUND ){
      /** if ( strcmp(instrume, "XIFU") ) { // Don't print these warnings for X-IFU event files.
	  SIXT_WARNING("'PHA' column not found! Falling back to 'signal' for spectra creation ...");
	  SIXT_WARNING("The spectrum will not be calibrated. ");
	  } **/
      fits_clear_errmsg();
      status = EXIT_SUCCESS;
      fits_get_colnum(ef, CASEINSEN, "signal", &csignal, &status);
      CHECK_STATUS_BREAK_WITH_FITSERROR(status);
      usesignal = 1;
    }

    CHECK_STATUS_BREAK_WITH_FITSERROR(status);

    if( usesignal==0 && par.usepha == 0 ){
      fits_get_colnum(ef, CASEINSEN, "PI", &coltmp, &status);
      // Print a warning if we find a PHA column in an X-IFU event file.
      if ( status == EXIT_SUCCESS && !strcmp(instrume, "XIFU") ) {
	SIXT_WARNING("Event file contains a PI column and is not compatible with this version of makespec");
      }

      if( status==COL_NOT_FOUND || status==KEY_NO_EXIST ){
	if ( strcmp(instrume, "XIFU") ) { // Don't print these warnings for X-IFU event files.
	  SIXT_WARNING("'PI' column not found! Falling back to 'PHA' for spectra creation ...");
	  SIXT_WARNING("The spectrum will not be calibrated. ");
	}
	// set default rmf file to nominal rmf used for simulation
	strcpy(rmffile, respfile );

	fits_clear_errmsg();
	status = EXIT_SUCCESS;
      } else {
	CHECK_STATUS_BREAK_WITH_FITSERROR(status);

	// now see if we find the PHA2PI information in the header used for pha2pi correction
	fits_read_key_longstr(ef, "PHA2PI", &pha2pi, NULL, &status); // ffgkls wants char ** longstr type

	if( status==VALUE_UNDEFINED || status==COL_NOT_FOUND || status==KEY_NO_EXIST){
	  headas_chat(5," *** warning : 'PHA2PI' key not found, but 'PI' column exits! Using 'PHA' values for spectra creation ...\n");
	  fits_clear_errmsg();
	  status = EXIT_SUCCESS;

	} else {
	  // now we have the PI column and the correction, so we use it
	  csignal = coltmp;

	  // set default rmf file to the Monte Carlo PI RMF if PIRMF key is available
	  if( strcasecmp("NONE",pirmf)!=0 ){
	    strcpy( rmffile, pirmf );
	  } else {
	    strcpy( rmffile, respfile );
	    headas_chat(3," *** warning : 'PIRMF' key not found, but analyzing 'PI' column! Using nominal RMF used for simulation instead ...\n");
	    headas_chat(3,"\t\t ... Model deviations expected!\n");
	  }
	}
	CHECK_STATUS_BREAK_WITH_FITSERROR(status);
      }
    }

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

    // we first check whether the user demands a different rmf or/and arf:
    // Check the rmf:
    char spcrespfile[MAXFILENAME];
    if (strcasecmp("NONE",par.RMFfile)){ // User demands a different rmf
      strcpy(spcrespfile,par.RMFfile);
      if( strcmp(spcrespfile,rmffile)!=0 ){
	headas_chat(5," *** WARNING: User RMF differs from expected RMF='%s'!\n",rmffile);
      }
    } else {                         // We use the same rmf as in the simulation
      strcpy(spcrespfile,rmffile);
    }

    // Decide which ARF to use: [USER] ARFfile, [fitshdr] ANCRfile, or SPECARF
    char arffile[MAXMSG];
    if( strcasecmp("NONE",specarf)!=0 ){
      strcpy( arffile, specarf );
    } else {
      strcpy( arffile, ancrfile );
    }

    // Check the user arf:
    char spcancrfile[MAXFILENAME];
    if (strcasecmp("NONE",par.ARFfile)){ // User demands a different arf
      strcpy(spcancrfile,par.ARFfile);
      if( strcmp(spcancrfile,arffile)!=0 ){
	headas_chat(5," *** WARNING: User ARF differs from expected ARF='%s'!\n",arffile);
      }
    } else {                         // We use the same arf as in the simulation
      strcpy(spcancrfile,arffile);
    }

    // We put the paths to the rmf and the arf into
    // resppathname and ancrpathname:
    char resppathname[2*MAXFILENAME];
    char ancrpathname[2*MAXFILENAME];
    if (strlen(par.RSPPath)>0) {
      // rmf:
      strcpy(resppathname, par.RSPPath);
      strcat(resppathname, "/");
      strcat(resppathname, spcrespfile);
      // arf:
      strcpy(ancrpathname, par.RSPPath);
      strcat(ancrpathname, "/");
      strcat(ancrpathname, spcancrfile);
    } else {
      // The file should be located in the working directory.
      strcpy(resppathname, spcrespfile);
      strcpy(ancrpathname, spcancrfile);
    }


    // Load the EBOUNDS of the RMF that will be used in the spectrum extraction.
    //   -> need to load full RMF to get FirstChannel correct
    //      (TODO: this can be simplified)
    struct RMF* rmf = loadRMF(resppathname,&status);
    CHECK_STATUS_BREAK(status);
    headas_chat(5," *** INFO : Using RMF='%s'!\n",resppathname);
    headas_chat(5," *** INFO : Using ARF='%s'!\n",ancrpathname);

    // If different rmf and/or arf are required, we need to check that the binning
    // is compatible with the ones used for the simulation:
    // Check the RMF (only necessary if using PHA/PI column, not if using signal
    // column):
    if (strcasecmp("NONE",par.RMFfile) && (usesignal != 1)){

      // Take the path to the rmf used in the simulation
      // we load the rmf used in the simulation:
      char simresppathname[2*MAXFILENAME];
      if (strlen(par.RSPPath)>0) {
      	strcpy(simresppathname, par.RSPPath);
      	strcat(simresppathname, "/");
      	strcat(simresppathname, respfile);
      } else {
      	// The file should be located in the working directory.
      	strcpy(simresppathname, respfile);
      }

      // Load the the RMF.
      struct RMF* simrmf=loadRMF(simresppathname,&status);
      CHECK_STATUS_BREAK(status);

      // We check the number of bins and the low and high energy:
      if((rmf->NumberChannels != simrmf->NumberChannels) ||
      	 (rmf->ChannelLowEnergy[0]   != simrmf->ChannelLowEnergy[0]) ||
      	 (rmf->ChannelHighEnergy[rmf->NumberChannels-2] != simrmf->ChannelHighEnergy[simrmf->NumberChannels-2])){
      	status=EXIT_FAILURE;
      	SIXT_ERROR("Required RMF has not the same binning as the original one");
        freeRMF(simrmf);
      	break;
      }

      // Check if FirstChannel is the same as in the RMF used for simulation.
      if (rmf->FirstChannel != simrmf->FirstChannel) {
        char msg[MAXFILENAME];
        snprintf(msg, sizeof(msg),
        		 "FirstChannel (=%li) of given RMF differs from FirstChannel (=%li) of RMF used for simulation",
        		 rmf->FirstChannel, simrmf->FirstChannel);
        status = EXIT_FAILURE;
        SIXT_ERROR(msg);
        freeRMF(simrmf);
        break;
      }
    }

    // Check the ARF:
    if (strcasecmp("NONE",par.ARFfile)){

      // Take the path to the arf used in the simulation
      // we load the arf used in the simulation:
      char simancrpathname[2*MAXFILENAME];
      if (strlen(par.RSPPath)>0) {
      	strcpy(simancrpathname, par.RSPPath);
      	strcat(simancrpathname, "/");
      	strcat(simancrpathname, ancrfile);
      } else {
      	// The file should be located in the working directory.
      	strcpy(simancrpathname, ancrfile);
      }

      // Load the ARFs.
      struct ARF* arf=loadARF(ancrpathname,&status);
      CHECK_STATUS_BREAK(status);
      struct ARF* simarf=loadARF(simancrpathname,&status);
      CHECK_STATUS_BREAK(status);

      // We check the number of bins and the low and high energy:
      if((arf->NumberEnergyBins != simarf->NumberEnergyBins) ||
      	 (arf->LowEnergy[0]     != simarf->LowEnergy[0]) ||
      	 (arf->HighEnergy[arf->NumberEnergyBins-2] != simarf->HighEnergy[simarf->NumberEnergyBins-2])){
      	status=EXIT_FAILURE;
      	SIXT_ERROR("Required ARF has not the same binning as the original one");
      	break;
      }
    }


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

    // Determine channel id from signal and rmf
    if( usesignal==1 ){
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
    }
    // Read channel id directly from PHA/PI
    else{
    	// LOOP over all events in the FITS table.
    	for (ii=0; ii<nrows; ii++) {

    		// Read the next event from the file.
    		long signal;
    		long fnull=0;
    		int anynul=0;
    		fits_read_col(ef, TLONG, csignal, ii+1, 1, 1,
    				&fnull, &signal, &anynul, &status);
    		CHECK_STATUS_BREAK(status);

    		// Add the event to the spectrum.
    		long idx=signal-rmf->FirstChannel;
    		if(idx>=0) {
    			assert(idx<rmf->NumberChannels);
    			spec[idx]++;
    		}
    	}
    	CHECK_STATUS_BREAK(status);
    	// END of loop over all events in the input file.
    }

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
    fits_update_key(sf, TSTRING, "LONGSTRN", "OGIP 1.0",
		    "The OGIP long string convention may be used", &status);
    fits_update_key_longstr(sf, "ANCRFILE", ancrpathname,
			    "ancillary response file", &status);
    fits_update_key_longstr(sf, "RESPFILE", resppathname,
			    "response file", &status);
    fits_update_key(sf, TSTRING, "BACKFILE", "",
		    "background file", &status);
    fits_update_key(sf, TLONG, "DETCHANS", &rmf->NumberChannels,
		    "number of detector channels", &status);
    if (pha2pi){
      fits_update_key_longstr(sf, "PHA2PI", pha2pi,
			      "Pha2Pi correction file", &status);
      free(pha2pi);
    }
    fits_update_key(sf, TSTRING, "CORRFILE", "NONE", "Correlation file", &status);
    fits_update_key(sf, TDOUBLE, "EXPOSURE", &exposure,
		    "exposure time", &status);
    fits_update_key(sf, TSTRING, "FilterExpr", par.EventFilter,
		    "used filter expression", &status);
    int logic=TRUE;
    fits_update_key(sf, TLOGICAL, "POISSERR", &logic, "Poissonian error", &status);
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

  status=ape_trad_query_file_name("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list file");
    return(status);
  }
  strcpy(par->EvtFile, sbuffer);
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

  //
  status=ape_trad_query_string("ANCRfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the ancilliary file");
    return(status);
  }
  strcpy(par->ARFfile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("RESPfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the response file");
    return(status);
  }
  strcpy(par->RMFfile, sbuffer);
  free(sbuffer);
  //

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
  status=ape_trad_query_int("usepha", &par->usepha);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the usepha parameter");
    return(status);
  }

  return(status);
}
