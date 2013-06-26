#include "runsixt.h"


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
  PhotonListFile* plf=NULL;

  // Impact list file.
  ImpactListFile* ilf=NULL;

  // Event list file.
  EventListFile* elf=NULL;

  // Pattern list file.
  PatternFile* patf=NULL;

  // Output file for progress status.
  FILE* progressfile=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("runsixt");
  set_toolversion("0.15");


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
    
    // Determine the event list output file.
    char eventlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.EventList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(eventlist_filename, par.Prefix);
      strcat(eventlist_filename, "events.fits");
    } else {
      strcpy(eventlist_filename, par.Prefix);
      strcat(eventlist_filename, par.EventList);
    }

    // Determine the pattern list output file.
    char patternlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.PatternList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(patternlist_filename, par.Prefix);
      strcat(patternlist_filename, "pattern.fits");
    } else {
      strcpy(patternlist_filename, par.Prefix);
      strcat(patternlist_filename, par.PatternList);
    }

    // Determine the random number seed.
    int seed;
    if (-1!=par.Seed) {
      seed=par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed=(int)time(NULL);
    }

    // Initialize the random number generator.
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
    inst=loadGenInst(xml_filename, &status);
    CHECK_STATUS_BREAK(status);
    
    // Set the usage of the detector background according to
    // the respective program parameter.
    setGenDetIgnoreBkg(inst->det, !par.Background);

    // Set up the Attitude.
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if ((strlen(par.Attitude)==0)||(0==strcmp(ucase_buffer, "NONE"))) {
      // Set up a simple pointing attitude.

      // First allocate memory.
      ac=getAttitude(&status);
      CHECK_STATUS_BREAK(status);

      ac->entry=(AttitudeEntry*)malloc(sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status = EXIT_FAILURE;
	SIXT_ERROR("memory allocation for Attitude failed");
	break;
      }

      // Set the values of the entries.
      ac->nentries=1;
      ac->entry[0]=defaultAttitudeEntry();
      ac->entry[0].time=par.TSTART;
      ac->entry[0].nz=unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);
      
      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      if ((ac->entry[0].time > par.TSTART) || 
	  (ac->entry[ac->nentries-1].time < par.TSTART+par.Exposure)) {
	status=EXIT_FAILURE;
	char msg[MAXMSG];
	sprintf(msg, "attitude data does not cover the "
		"specified period from %lf to %lf!", 
		par.TSTART, par.TSTART+par.Exposure);
	SIXT_ERROR(msg);
	break;
      }
    }
    // END of setting up the attitude.

    // Optional GTI file.
    if (strlen(par.GTIfile)>0) {
      strcpy(ucase_buffer, par.GTIfile);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	gti=loadGTI(par.GTIfile, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

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

    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    double tstop;
    if (NULL==gti) {
      tstop=par.TSTART+par.Exposure;
    } else {
      tstop=gti->stop[gti->nentries-1];
    }

    // Open the output photon list file.
    if (strlen(photonlist_filename)>0) {
      plf=openNewPhotonListFile(photonlist_filename, 
				telescop, instrume, "Normal",
				par.MJDREF, 0.0, par.TSTART, tstop,
				par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list file.
    if (strlen(impactlist_filename)>0) {
      ilf=openNewImpactListFile(impactlist_filename, 
				telescop, instrume, "Normal",
				par.MJDREF, 0.0, par.TSTART, tstop,
				par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list file.
    elf=openNewEventListFile(eventlist_filename, 
			     telescop, instrume, "Normal",
			     par.MJDREF, 0.0, par.TSTART, tstop,
			     inst->det->pixgrid->xwidth,
			     inst->det->pixgrid->ywidth,
			     par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Define the event list file as output file.
    setGenDetEventListFile(inst->det, elf);

    // Open the output pattern list file.
    patf=openNewPatternFile(patternlist_filename, 
			    telescop, instrume, "Normal",			    
			    par.MJDREF, 0.0, par.TSTART, tstop,
			    inst->det->pixgrid->xwidth,
			    inst->det->pixgrid->ywidth,
			    par.clobber, &status);
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
      ra *= 180./M_PI;
      dec*= 180./M_PI;
      rollangle*= 180./M_PI;

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
      if (NULL!=elf) {
	fits_update_key(elf->fptr, TDOUBLE, "RA_PNT", &ra,
			"RA of pointing direction [deg]", &status);
	fits_update_key(elf->fptr, TDOUBLE, "DEC_PNT", &dec,
			"Dec of pointing direction [deg]", &status);
	fits_update_key(elf->fptr, TFLOAT, "PA_PNT", &rollangle,
			"Roll angle [deg]", &status);
	CHECK_STATUS_BREAK(status);
      }

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
      if (NULL!=elf) {
	fits_update_key(elf->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
      }
      fits_update_key(patf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		      "attitude file", &status);
      CHECK_STATUS_BREAK(status);
    }

    // TLMIN and TLMAX of PI column.
    char keystr[MAXMSG];
    long value;
    if (NULL!=elf) {
      sprintf(keystr, "TLMIN%d", elf->cpi);
      value=inst->det->rmf->FirstChannel;
      fits_update_key(elf->fptr, TLONG, keystr, &value, "", &status);
      sprintf(keystr, "TLMAX%d", elf->cpi);
      value=inst->det->rmf->FirstChannel+inst->det->rmf->NumberChannels-1;
      fits_update_key(elf->fptr, TLONG, keystr, &value, "", &status);
      CHECK_STATUS_BREAK(status);
    }
    sprintf(keystr, "TLMIN%d", patf->cpi);
    value=inst->det->rmf->FirstChannel;
    fits_update_key(patf->fptr, TLONG, keystr, &value, "", &status);
    sprintf(keystr, "TLMAX%d", patf->cpi);
    value=inst->det->rmf->FirstChannel+inst->det->rmf->NumberChannels-1;
    fits_update_key(patf->fptr, TLONG, keystr, &value, "", &status);
    CHECK_STATUS_BREAK(status);

    // --- End of opening files ---


    // --- Simulation Process ---

    headas_chat(3, "start simulation ...\n");

    // Simulation progress status (running from 0 to 100).
    unsigned int progress=0;
    if (NULL==progressfile) {
      headas_chat(2, "\r%.1lf %%", 0.);
      fflush(NULL);
    } else {
      rewind(progressfile);
      fprintf(progressfile, "%.2lf", 0.);
      fflush(progressfile);	
    }

    // Determine the total length of the time interval to
    // be simulated.
    double totalsimtime=0.;
    double simtime=0.;
    if (NULL==gti) {
      totalsimtime=par.Exposure;
    } else {
      unsigned long ii; 
      for (ii=0; ii<gti->nentries; ii++) {
	totalsimtime+=gti->stop[ii]-gti->start[ii];
      }
    }


    // Current bin in the GTI collection.
    unsigned long gtibin=0;
    // Loop over all intervals in the GTI collection.
    do {
      // Currently regarded interval.
      double t0, t1;
      
      // Determine the currently regarded interval.
      if (NULL==gti) {
	t0=par.TSTART;
	t1=par.Exposure;
      } else {
	t0=gti->start[gtibin];
	t1=gti->stop[gtibin];
      }

      // Set the start time for the instrument model.
      setGenDetStartTime(inst->det, t0);

      // Loop over photon generation and processing
      // till the time of the photon exceeds the requested
      // time interval.
      do {

	// Photon generation.
	Photon ph;
	int isph=phgen(ac, srccat, MAX_N_SIMPUT, t0, t1,
		       par.MJDREF, par.dt, 
		       inst->tel->fov_diameter, &ph, &status);
	CHECK_STATUS_BREAK(status);

	// If no photon has been generated, break the loop.
	if (0==isph) break;

	// Check if the photon still is within the requested
	// exposre time.
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
	    headas_chat(2, "\r%.1lf %%", progress*1.);
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
      if (NULL!=gti) {
	simtime+=gti->stop[gtibin]-gti->start[gtibin];
	gtibin++;
	if (gtibin>=gti->nentries) break;
      }

    } while (NULL!=gti);
    CHECK_STATUS_BREAK(status);
    // End of loop over the individual GTI intervals.
    
      
    // Progress output.
    if (NULL==progressfile) {
      headas_chat(2, "\r%.1lf %%\n", 100.);
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
      copyEvents2PatternFile(elf, patf, &status);
      CHECK_STATUS_BREAK(status);
    }
    
    // Close files in order to save memory.
    freePhotonListFile(&plf, &status);
    freeImpactListFile(&ilf, &status);
    freeEventListFile(&elf, &status);

    // Run the event projection.
    headas_chat(3, "start sky projection ...\n");
    phproj(inst, ac, patf, par.TSTART, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  destroyPatternFile(&patf, &status);
  freeEventListFile(&elf, &status);
  freeImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);
  for (ii=0; ii<MAX_N_SIMPUT; ii++) {
    freeSourceCatalog(&(srccat[ii]), &status);
  }
  freeGTI(&gti);
  freeAttitude(&ac);
  destroyGenInst(&inst, &status);

  if (NULL!=progressfile) {
    fclose(progressfile);
    progressfile=NULL;
  }

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
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

  status=ape_trad_query_string("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pattern list");
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
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

  status=ape_trad_query_bool("Background", &par->Background);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the background flag");
    return(status);
  }

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the attitude");
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the right ascension of the telescope "
	       "pointing");
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the declination of the telescope "
	       "pointing");
    return(status);
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

  status=ape_trad_query_int("seed", &par->Seed);
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


