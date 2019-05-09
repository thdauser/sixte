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

#include "runsixt.h"
#include "parinput.h"


int runsixt_main()
{
  // Program parameters.
  struct Parameters par;

  // Instrument setup.
  GenInst* inst=NULL;

  // Attitude.
  Attitude* ac=NULL;

  // GTI collection.
  GTI* gti=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat[MAX_N_SIMPUT];
  unsigned long ii;
  for (ii=0; ii<MAX_N_SIMPUT; ii++) {
    srccat[ii]=NULL;
  }

  // Photon list file.
  PhotonFile* plf=NULL;

  // Impact list file.
  ImpactFile* ilf=NULL;

  // Single-pixel event file.
  EventFile* elf=NULL;

  // Pattern event file.
  EventFile* patf=NULL;

  // Output file for progress status.
  FILE* progressfile=NULL;

  // Pha2Pi correction file
  Pha2Pi* p2p=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("runsixt");
  set_toolversion("0.19");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=runsixt_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Determine the prefix for the output files.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Prefix);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(par.Prefix, "");
    }

    // Determine the photon list output file.
    char photonlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.PhotonList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(photonlist_filename, "");
    } else {
      strcpy(photonlist_filename, par.Prefix);
      strcat(photonlist_filename, par.PhotonList);
    }

    // Determine the impact list output file.
    char impactlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.ImpactList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(impactlist_filename, "");
    } else {
      strcpy(impactlist_filename, par.Prefix);
      strcat(impactlist_filename, par.ImpactList);
    }

    // Determine the single-pixel event output file.
    char rawdata_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.RawData);
    strtoupper(ucase_buffer);
    int delete_rawdata = 0;
    if (0==strcmp(ucase_buffer,"NONE")) {
      delete_rawdata = 1;
      strcpy(rawdata_filename, par.Prefix);
      strcat(rawdata_filename, "raw.fits");
    } else {
      strcpy(rawdata_filename, par.Prefix);
      strcat(rawdata_filename, par.RawData);
    }

    // Determine the event pattern output file.
    char evtfile_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.EvtFile);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(evtfile_filename, par.Prefix);
      strcat(evtfile_filename, "evt.fits");
    } else {
      strcpy(evtfile_filename, par.Prefix);
      strcat(evtfile_filename, par.EvtFile);
    }

    // Initialize the random number generator.
    unsigned int seed=getSeed(par.Seed);
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Set the progress status output file.
    strcpy(ucase_buffer, par.ProgressFile);
    strtoupper(ucase_buffer);
    if (0!=strcmp(ucase_buffer, "STDOUT")) {
      progressfile=fopen(par.ProgressFile, "w+");
      char msg[MAXMSG];
      sprintf(msg, "could not open file '%s' for output of progress status",
	      par.ProgressFile);
      CHECK_NULL_BREAK(progressfile, status, msg);
    }

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile,
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, seed, &status);
    CHECK_STATUS_BREAK(status);

    // Initialize & load Pha2Pi File (NULL if not set)
    p2p = initPha2Pi_from_GenInst( inst, seed, &status);
    CHECK_STATUS_BREAK_WITH_FITSERROR(status);

    // Set the usage of the detector background according to
    // the respective program parameter.
    setGenDetIgnoreBkg(inst->det, !par.Background);

    // Set up the Attitude.
    if (par.Attitude=='\0' || (strlen(par.Attitude)==0) ) {
      // Set up a pointing attitude.
      ac=getPointingAttitude(par.MJDREF, par.TSTART, par.TSTART+par.Exposure,
			     par.RA*M_PI/180., par.Dec*M_PI/180., par.rollangle*M_PI/180., &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval for the simulation
      // is a subset of the period described by the attitude file.
      checkAttitudeTimeCoverage(ac, par.MJDREF, par.TSTART,
				par.TSTART+par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of setting up the attitude.

    // Get a GTI.
    gti=getGTIFromFileOrContinuous(par.GTIfile,
				   par.TSTART, par.TSTART+par.Exposure,
				   par.MJDREF, &status);
    CHECK_STATUS_BREAK(status);

    // Load the SIMPUT X-ray source catalogs.
    srccat[0]=loadSourceCatalog(par.Simput, inst->tel->arf, &status);
    CHECK_STATUS_BREAK(status);

    // Optional 2nd catalog.
    if (strlen(par.Simput2)>0) {
      strcpy(ucase_buffer, par.Simput2);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[1]=loadSourceCatalog(par.Simput2, inst->tel->arf, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 3rd catalog.
    if (strlen(par.Simput3)>0) {
      strcpy(ucase_buffer, par.Simput3);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[2]=loadSourceCatalog(par.Simput3, inst->tel->arf, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 4th catalog.
    if (strlen(par.Simput4)>0) {
      strcpy(ucase_buffer, par.Simput4);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[3]=loadSourceCatalog(par.Simput4, inst->tel->arf, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 5th catalog.
    if (strlen(par.Simput5)>0) {
      strcpy(ucase_buffer, par.Simput5);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[4]=loadSourceCatalog(par.Simput5, inst->tel->arf, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 6th catalog.
    if (strlen(par.Simput6)>0) {
      strcpy(ucase_buffer, par.Simput6);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[5]=loadSourceCatalog(par.Simput6, inst->tel->arf, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // --- End of Initialization ---


    // --- Open and set up files ---

    char telescop[MAXMSG]={""}, instrume[MAXMSG]={""}, filter[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    if (NULL!=inst->tel->arf->Filter) {
      strcpy(filter, inst->tel->arf->Filter);
    }
    double tstop=gti->stop[gti->ngti-1];

    // Open the output photon list file.
    if (strlen(photonlist_filename)>0) {
      plf=openNewPhotonFile(photonlist_filename,
			    telescop, instrume, filter,
			    inst->tel->arf_filename, inst->det->rmf_filename,
			    par.MJDREF, 0.0, par.TSTART, tstop,
			    par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list file.
    if (strlen(impactlist_filename)>0) {
      ilf=openNewImpactFile(impactlist_filename,
			    telescop, instrume, filter,
			    inst->tel->arf_filename, inst->det->rmf_filename,
			    par.MJDREF, 0.0, par.TSTART, tstop,
			    par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list file.
    elf=openNewEventFile(rawdata_filename,
			 telescop, instrume, filter,
			 inst->tel->arf_filename, inst->det->rmf_filename,
			 par.MJDREF, 0.0, par.TSTART, tstop,
			 inst->det->pixgrid->xwidth,
			 inst->det->pixgrid->ywidth,
			 par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Define the event file as output file.
    setGenDetEventFile(inst->det, elf);

    // Open the output pattern list file.
    patf=openNewEventFile(evtfile_filename,
			  telescop, instrume, filter,
			  inst->tel->arf_filename, inst->det->rmf_filename,
			  par.MJDREF, 0.0, par.TSTART, tstop,
			  inst->det->pixgrid->xwidth,
			  inst->det->pixgrid->ywidth,
			  par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    float rotation_angle=inst->det->pixgrid->rota*180./M_PI;
    fits_update_key(elf->fptr, TFLOAT, "CCDROTA", &rotation_angle, "CCD rotation angle [deg]", &status);
    fits_update_key(patf->fptr, TFLOAT, "CCDROTA", &rotation_angle, "CCD rotation angle [deg]", &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    // If this is a pointing attitude, store the direction in the output
    // photon list.
    if (1==ac->nentries) {
      // Determine the telescope pointing direction and roll angle.
      Vector pointing=getTelescopeNz(ac, par.TSTART, &status);
      CHECK_STATUS_BREAK(status);

      // Direction.
      double ra, dec;
      calculate_ra_dec(pointing, &ra, &dec);

      // Roll angle.
      float rollangle=getRollAngle(ac, par.TSTART, &status);
      CHECK_STATUS_BREAK(status);

      // Store the RA and Dec information in the FITS header.
      ra *=180./M_PI;
      dec*=180./M_PI;
      rollangle*=180./M_PI;

      // Photon list file.
      if (NULL!=plf) {
	fits_update_key(plf->fptr, TDOUBLE, "RA_PNT", &ra,
			"RA of pointing direction [deg]", &status);
	fits_update_key(plf->fptr, TDOUBLE, "DEC_PNT", &dec,
			"Dec of pointing direction [deg]", &status);
	fits_update_key(plf->fptr, TFLOAT, "PA_PNT", &rollangle,
			"Roll angle [deg]", &status);
	CHECK_STATUS_BREAK(status);
      }

      // Impact list file.
      if (NULL!=ilf) {
	fits_update_key(ilf->fptr, TDOUBLE, "RA_PNT", &ra,
			"RA of pointing direction [deg]", &status);
	fits_update_key(ilf->fptr, TDOUBLE, "DEC_PNT", &dec,
			"Dec of pointing direction [deg]", &status);
	fits_update_key(ilf->fptr, TFLOAT, "PA_PNT", &rollangle,
			"Roll angle [deg]", &status);
	CHECK_STATUS_BREAK(status);
      }

      // Event list file.
      fits_update_key(elf->fptr, TDOUBLE, "RA_PNT", &ra,
		      "RA of pointing direction [deg]", &status);
      fits_update_key(elf->fptr, TDOUBLE, "DEC_PNT", &dec,
		      "Dec of pointing direction [deg]", &status);
      fits_update_key(elf->fptr, TFLOAT, "PA_PNT", &rollangle,
		      "Roll angle [deg]", &status);
      CHECK_STATUS_BREAK(status);

      // Pattern list file.
      fits_update_key(patf->fptr, TDOUBLE, "RA_PNT", &ra,
		      "RA of pointing direction [deg]", &status);
      fits_update_key(patf->fptr, TDOUBLE, "DEC_PNT", &dec,
		      "Dec of pointing direction [deg]", &status);
      fits_update_key(patf->fptr, TFLOAT, "PA_PNT", &rollangle,
		      "Roll angle [deg]", &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // An explicit attitude file is given.
      if (NULL!=plf) {
	fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
      }
      if (NULL!=ilf) {
	fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
      }
      fits_update_key(elf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		      "attitude file", &status);
      fits_update_key(patf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		      "attitude file", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Event type.
    fits_update_key(elf->fptr, TSTRING, "EVTYPE", "PIXEL",
		    "event type", &status);
    CHECK_STATUS_BREAK(status);

    // TLMIN and TLMAX of PI column.
    char keystr[MAXMSG];
    long value;
    sprintf(keystr, "TLMIN%d", elf->cpha);
    value=inst->det->rmf->FirstChannel;
    fits_update_key(elf->fptr, TLONG, keystr, &value, "", &status);
    sprintf(keystr, "TLMAX%d", elf->cpha);
    value=inst->det->rmf->FirstChannel+inst->det->rmf->NumberChannels-1;
    fits_update_key(elf->fptr, TLONG, keystr, &value, "", &status);
    CHECK_STATUS_BREAK(status);

    sprintf(keystr, "TLMIN%d", patf->cpha);
    value=inst->det->rmf->FirstChannel;
    fits_update_key(patf->fptr, TLONG, keystr, &value, "", &status);
    sprintf(keystr, "TLMAX%d", patf->cpha);
    value=inst->det->rmf->FirstChannel+inst->det->rmf->NumberChannels-1;
    fits_update_key(patf->fptr, TLONG, keystr, &value, "", &status);
    CHECK_STATUS_BREAK(status);

    // --- End of opening files ---


    // --- Simulation Process ---

    headas_chat(3, "start simulation ...\n");

    // Simulation progress status (running from 0 to 100).
    unsigned int progress=0;
    if (NULL==progressfile) {
      headas_chat(2, "\r%.0lf %%", 0.);
      fflush(NULL);
    } else {
      rewind(progressfile);
      fprintf(progressfile, "%.2lf", 0.);
      fflush(progressfile);
    }

    // Determine the total length of the time interval to
    // be simulated.
    double totalsimtime=sumGTI(gti);
    fits_update_key(patf->fptr, TDOUBLE, "EXPOSURE", &totalsimtime,
		    "exposure time [s]", &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all intervals in the GTI collection.
    int gtibin=0;
    double simtime=0.;
    do {
    	// Currently regarded interval.
    	double t0=gti->start[gtibin];
    	double t1=gti->stop[gtibin];

    	// Set the start time for the instrument model.
    	setGenDetStartTime(inst->det, t0);

    	// Loop over photon generation and processing
    	// till the time of the photon exceeds the requested
    	// time interval.
    	do {

    		// Photon generation.
    		Photon ph;
    		int isph=phgen(ac, srccat, MAX_N_SIMPUT, t0, t1, par.MJDREF, par.dt,
    				inst->tel->fov_diameter, &ph, &status);
    		CHECK_STATUS_BREAK(status);

    		// If no photon has been generated, break the loop.
    		if (0==isph) break;

    		// Check if the photon still is within the requested
    		// exposure time.
    		assert(ph.time<=t1);

    		// If requested, write the photon to the output file.
    		if (NULL!=plf) {
    			status=addPhoton2File(plf, &ph);
    			CHECK_STATUS_BREAK(status);
    		}

    		// Photon imaging.
    		Impact imp;
    		int isimg=phimg(inst->tel, ac, &ph, &imp, &status);
    		CHECK_STATUS_BREAK(status);

    		// If the photon is not imaged but lost in the optical system,
    		// continue with the next one.
    		if (0==isimg) continue;

    		// If requested, write the impact to the output file.
    		if (NULL!=ilf) {
    			addImpact2File(ilf, &imp, &status);
    			CHECK_STATUS_BREAK(status);
    		}

    		// Photon Detection.
    		phdetGenDet(inst->det, &imp, t1, &status);
    		CHECK_STATUS_BREAK(status);

    		// Program progress output.
    		while((unsigned int)((ph.time-t0+simtime)*100./totalsimtime)>progress) {
    			progress++;
    			if (NULL==progressfile) {
    				headas_chat(2, "\r%.0lf %%", progress*1.);
    				fflush(NULL);
    			} else {
    				rewind(progressfile);
    				fprintf(progressfile, "%.2lf", progress*1./100.);
    				fflush(progressfile);
    			}
    		}

    	} while(1);
    	CHECK_STATUS_BREAK(status);
    	// END of photon processing loop for the current interval.

    	// Clear the detector.
    	phdetGenDet(inst->det, NULL, t1, &status);
    	CHECK_STATUS_BREAK(status);
    	long jj;
    	for(jj=0; jj<inst->det->pixgrid->ywidth; jj++) {
    		GenDetClearLine(inst->det, jj);
    	}

    	// Proceed to the next GTI interval.
    	simtime+=gti->stop[gtibin]-gti->start[gtibin];
    	gtibin++;
    	if (gtibin>=gti->ngti) break;

    } while (1);
    CHECK_STATUS_BREAK(status);
    // End of loop over the individual GTI intervals.


    // Progress output.
    if (NULL==progressfile) {
    	headas_chat(2, "\r%.0lf %%\n", 100.);
    	fflush(NULL);
    } else {
    	rewind(progressfile);
    	fprintf(progressfile, "%.2lf", 1.);
    	fflush(progressfile);
    }

    // Perform a pattern analysis, only if split events are simulated.
    if (GS_NONE!=inst->det->split->type) {
    	// Pattern analysis.
    	headas_chat(3, "start event pattern analysis ...\n");
    	phpat(inst->det, elf, patf, par.SkipInvalids, &status);
    	CHECK_STATUS_BREAK(status);

    } else {
    	// If no split events are simulated, simply copy the event list
    	// to a pattern list.
    	headas_chat(3, "copy events to pattern file ...\n");
    	copyEventFile(elf, patf,
    			inst->det->threshold_event_lo_keV,
				inst->det->threshold_pattern_up_keV,
				&status);
    	CHECK_STATUS_BREAK(status);
    	fits_update_key(patf->fptr, TSTRING, "EVTYPE", "PATTERN",
    			"event type", &status);
    	CHECK_STATUS_BREAK(status);
    }

    // Store the GTI extension in the event file.
    saveGTIExt(elf->fptr, "STDGTI", gti, &status);
    CHECK_STATUS_BREAK(status);

    // Close files in order to save memory.
    freePhotonFile(&plf, &status);
    freeImpactFile(&ilf, &status);
    freeEventFile(&elf, &status);

    // Run the event projection.
    headas_chat(3, "start sky projection ...\n");
    phproj(inst, ac, patf, par.TSTART, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

    // Store the GTI extension in the pattern file.
    saveGTIExt(patf->fptr, "STDGTI", gti, &status);
    CHECK_STATUS_BREAK(status);

    // Run PI correction on Pattern file.
    if( p2p != NULL ){
    	headas_chat(3, "start Pha2Pi correction ...\n");
    	pha2pi_correct_eventfile( patf, p2p, inst->filepath, inst->det->rmf_filename, &status);
    	CHECK_STATUS_BREAK_WITH_FITSERROR(status);
    }

    // --- End of simulation process ---
    // remove RawData files if not requested
    if (delete_rawdata){
    	headas_chat(5,"removing unwanted RawData file %s \n",rawdata_filename);
    	status = remove (rawdata_filename);
    	CHECK_STATUS_BREAK(status);
    }


  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---

  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeEventFile(&patf, &status);
  freeEventFile(&elf, &status);
  freeImpactFile(&ilf, &status);
  freePhotonFile(&plf, &status);
  for (ii=0; ii<MAX_N_SIMPUT; ii++) {
    freeSourceCatalog(&(srccat[ii]), &status);
  }
  freePha2Pi(&p2p);
  freeGTI(&gti);
  freeAttitude(&ac);
  destroyGenInst(&inst, &status);

  if (NULL!=progressfile) {
    fclose(progressfile);
    progressfile=NULL;
  }

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int runsixt_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("Prefix", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the prefix for the output files");
    return(status);
  }
  strcpy(par->Prefix, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("PhotonList", &sbuffer);
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

  status=ape_trad_query_string("RawData", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list");
    return(status);
  }
  strcpy(par->RawData, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pattern list");
    return(status);
  }
  strcpy(par->EvtFile, sbuffer);
  free(sbuffer);

  query_simput_parameter_file_name("XMLFile", &(par->XMLFile), &status);

  // only load Mission, Instrument and Mode if XMLFile is not given
  if (par->XMLFile=='\0'){
    query_simput_parameter_string("Mission", &(par->Mission), &status);
    query_simput_parameter_string("Instrument", &(par->Instrument), &status);
    query_simput_parameter_string("Mode", &(par->Mode), &status);
  } else {
	  // set to default values
    par->Mission = NULL;
    par->Instrument = NULL;
    par->Mode = NULL;
	  headas_chat(5, "using xml file: %s \n", par->XMLFile);
  }

  status=ape_trad_query_bool("Background", &par->Background);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the background flag");
    return(status);
  }

  query_simput_parameter_file_name("Attitude", &(par->Attitude), &status);

  // only load RA,Dec if Attitude is not given
  if (par->Attitude=='\0'){
	  query_simput_parameter_float("RA",&(par->RA),&status);
	  query_simput_parameter_float("Dec",&(par->Dec),&status);
	  query_simput_parameter_float("rollangle",&(par->rollangle),&status);
  } else {
	  // set to default values
	  par->RA=0.0;
	  par->Dec=0.0;
	  par->rollangle=0.0;
	  headas_chat(3, "using Attiude File: %s \n",par->Attitude);
  }

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the SIMPUT file");
    return(status);
  }
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Simput2", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the second SIMPUT file");
    return(status);
  }
  strcpy(par->Simput2, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Simput3", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the third SIMPUT file");
    return(status);
  }
  strcpy(par->Simput3, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Simput4", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the forth SIMPUT file");
    return(status);
  }
  strcpy(par->Simput4, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Simput5", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the fifth SIMPUT file");
    return(status);
  }
  strcpy(par->Simput5, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Simput6", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the sixth SIMPUT file");
    return(status);
  }
  strcpy(par->Simput6, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("GTIfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the GTI file");
    return(status);
  }
  strcpy(par->GTIfile, sbuffer);
  free(sbuffer);

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

  status=ape_trad_query_double("dt", &par->dt);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading dt");
    return(status);
  }

  status=ape_trad_query_bool("SkipInvalids", &par->SkipInvalids);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the SkipInvalids parameter");
    return(status);
  }

  status=ape_trad_query_int("Seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the seed for the random number generator");
    return(status);
  }

  status=ape_trad_query_string("ProgressFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the progress status file");
    return(status);
  }
  strcpy(par->ProgressFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}
