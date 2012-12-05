#include "erosim.h"


int erosim_main() 
{
  // Program parameters.
  struct Parameters par;
  
  // Individual sub-instruments.
  GenInst* subinst[7]={NULL, NULL, NULL, NULL, NULL, NULL, NULL};

  // Fake telescope ARF with the composite effective area of all
  // 7 sub-telescopes. This is used for the photon generation
  // from the SIMPUT catalog.
  struct ARF* arf7=NULL;
  
  // FoV of the combined 7 telescopes. If there is not misalignment,
  // it corresponds to the FoV of a single sub-telescope.
  float fov7;

  // Attitude.
  AttitudeCatalog* ac=NULL;

  // GTI collection.
  GTI* gti=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat[MAX_N_SIMPUT];
  unsigned long ii;
  for (ii=0; ii<MAX_N_SIMPUT; ii++) {
    srccat[ii]=NULL;
  }

  // Photon list files.
  PhotonListFile* plf[7]={NULL, NULL, NULL, NULL, NULL, NULL, NULL};

  // Impact list file.
  ImpactListFile* ilf[7]={NULL, NULL, NULL, NULL, NULL, NULL, NULL};

  // Event list file.
  EventListFile* elf[7]={NULL, NULL, NULL, NULL, NULL, NULL, NULL};

  // Pattern list file.
  PatternFile* patf[7]={NULL, NULL, NULL, NULL, NULL, NULL, NULL};

  // Output file for progress status.
  FILE* progressfile=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("erosim");
  set_toolversion("0.02");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=erosim_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Determine the prefix for the output files.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Prefix);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(par.Prefix, "");
    }

    // Determine the photon list output files.
    char photonlist_filename_template[MAXFILENAME];
    strcpy(ucase_buffer, par.PhotonList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(photonlist_filename_template, "");
    } else {
      strcpy(photonlist_filename_template, par.Prefix);
      strcat(photonlist_filename_template, "tel%d_");
      strcat(photonlist_filename_template, par.PhotonList);
    }

    // Determine the impact list output file.
    char impactlist_filename_template[MAXFILENAME];
    strcpy(ucase_buffer, par.ImpactList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(impactlist_filename_template, "");
    } else {
      strcpy(impactlist_filename_template, par.Prefix);
      strcat(impactlist_filename_template, "tel%d_");
      strcat(impactlist_filename_template, par.ImpactList);
    }
    
    // Determine the event list output file.
    char eventlist_filename_template[MAXFILENAME];
    strcpy(ucase_buffer, par.EventList);
    strtoupper(ucase_buffer);
    strcpy(eventlist_filename_template, par.Prefix);
    strcat(eventlist_filename_template, "ccd%d_");
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcat(eventlist_filename_template, "events.fits");
    } else {
      strcat(eventlist_filename_template, par.EventList);
    }

    // Determine the pattern list output file.
    char patternlist_filename_template[MAXFILENAME];
    strcpy(ucase_buffer, par.PatternList);
    strtoupper(ucase_buffer);
    strcpy(patternlist_filename_template, par.Prefix);
    strcat(patternlist_filename_template, "ccd%d_");
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcat(patternlist_filename_template, "pattern.fits");
    } else {
      strcat(patternlist_filename_template, par.PatternList);
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
		     "SRG", "eROSITA", "", &status);
    CHECK_STATUS_BREAK(status);

    // Load the configurations of all seven sub-instruments.
    for (ii=0; ii<7; ii++) {
      // Check if a particular XML file is given for this 
      // sub-instrument. If not, use the default XML file.
      char buffer[MAXFILENAME];
      char ubuffer[MAXFILENAME];
      switch (ii) {
      case 0:
	strcpy(buffer, par.XMLFile1);
	break;
      case 1:
	strcpy(buffer, par.XMLFile2);
	break;
      case 2:
	strcpy(buffer, par.XMLFile3);
	break;
      case 3:
	strcpy(buffer, par.XMLFile4);
	break;
      case 4:
	strcpy(buffer, par.XMLFile5);
	break;
      case 5:
	strcpy(buffer, par.XMLFile6);
	break;
      case 6:
	strcpy(buffer, par.XMLFile7);
	break;
      default:
	break;
      }
      strcpy(ubuffer, buffer);
      strtoupper(ubuffer);
      if (0==strcmp(ubuffer, "NONE")) {
	strcpy(buffer, xml_filename);
      }

      // Load the instrument configuration either with the
      // specific (if available) or the default XML file.
      subinst[ii]=loadGenInst(buffer, &status);
      CHECK_STATUS_BREAK(status);

      // Set the usage of the detector background according to
      // the respective program parameter.
      setGenDetIgnoreBkg(subinst[ii]->det, !par.Background);
    }
    CHECK_STATUS_BREAK(status);
    
    // Determine a fake ARF with the combined effective area of
    // all seven sub-telescopes.
    arf7=(struct ARF*)malloc(sizeof(struct ARF));
    CHECK_NULL_BREAK(arf7, status, "memory allocation for fake ARF failed");
    arf7->NumberEnergyBins=subinst[0]->tel->arf->NumberEnergyBins;
    arf7->LowEnergy= (float*)malloc(arf7->NumberEnergyBins*sizeof(float));
    arf7->HighEnergy=(float*)malloc(arf7->NumberEnergyBins*sizeof(float));
    arf7->EffArea=   (float*)malloc(arf7->NumberEnergyBins*sizeof(float));
    long kk;
    for (kk=0; kk<arf7->NumberEnergyBins; kk++) {
      arf7->LowEnergy[kk] =subinst[0]->tel->arf->LowEnergy[kk];
      arf7->HighEnergy[kk]=subinst[0]->tel->arf->HighEnergy[kk];
      arf7->EffArea[kk]   =subinst[0]->tel->arf->EffArea[kk] *7.; // Factor 7 !
    }
    strcpy(arf7->ARFVersion, subinst[0]->tel->arf->ARFVersion);
    strcpy(arf7->Telescope , subinst[0]->tel->arf->Telescope );
    strcpy(arf7->Instrument, subinst[0]->tel->arf->Instrument);
    strcpy(arf7->Detector  , subinst[0]->tel->arf->Detector  );
    strcpy(arf7->Filter    , subinst[0]->tel->arf->Filter    );
    strcpy(arf7->ARFExtensionName, subinst[0]->tel->arf->ARFExtensionName);
    
    // The FoV is the same as for an individual sub-telescope.
    // We ignore any misalignment for the moment.
    fov7=subinst[0]->tel->fov_diameter;
    
    // Set up the Attitude.
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if ((strlen(par.Attitude)==0)||(0==strcmp(ucase_buffer, "NONE"))) {
      // Set up a simple pointing attitude.

      // First allocate memory.
      ac=getAttitudeCatalog(&status);
      CHECK_STATUS_BREAK(status);

      ac->entry=(AttitudeEntry*)malloc(sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status = EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeCatalog failed");
	break;
      }

      // Set the values of the entries.
      ac->nentries=1;
      ac->entry[0]=defaultAttitudeEntry();
      ac->entry[0].time=par.TSTART;
      ac->entry[0].nz=unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
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
    srccat[0]=loadSourceCatalog(par.Simput, arf7, &status);
    CHECK_STATUS_BREAK(status);

    // Optional 2nd catalog.
    if (strlen(par.Simput2)>0) {
      strcpy(ucase_buffer, par.Simput2);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[1]=loadSourceCatalog(par.Simput2, arf7, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 3rd catalog.
    if (strlen(par.Simput3)>0) {
      strcpy(ucase_buffer, par.Simput3);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[2]=loadSourceCatalog(par.Simput3, arf7, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 4th catalog.
    if (strlen(par.Simput4)>0) {
      strcpy(ucase_buffer, par.Simput4);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[3]=loadSourceCatalog(par.Simput4, arf7, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 5th catalog.
    if (strlen(par.Simput5)>0) {
      strcpy(ucase_buffer, par.Simput5);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[4]=loadSourceCatalog(par.Simput5, arf7, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 6th catalog.
    if (strlen(par.Simput6)>0) {
      strcpy(ucase_buffer, par.Simput6);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[5]=loadSourceCatalog(par.Simput6, arf7, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // --- End of Initialization ---


    // --- Open and set up files ---

    // Open the output photon list files.
    if (strlen(photonlist_filename_template)>0) {
      for (ii=0; ii<7; ii++) {
	char photonlist_filename[MAXFILENAME];
	sprintf(photonlist_filename, photonlist_filename_template, ii);
	plf[ii]=openNewPhotonListFile(photonlist_filename, par.clobber, &status);
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list files.
    if (strlen(impactlist_filename_template)>0) {
      for (ii=0; ii<7; ii++) {
	char impactlist_filename[MAXFILENAME];
	sprintf(impactlist_filename, impactlist_filename_template, ii);
	ilf[ii]=openNewImpactListFile(impactlist_filename, par.clobber, &status);
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list files.
    for (ii=0; ii<7; ii++) {
      // Open the file.
      char eventlist_filename[MAXFILENAME];
      sprintf(eventlist_filename, eventlist_filename_template, ii);
      elf[ii]=openNewEventListFile(eventlist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);

      // Define the event list file as output file for the respective
      // detector.
      setGenDetEventListFile(subinst[ii]->det, elf[ii]);
    }
    CHECK_STATUS_BREAK(status);

    // Open the output pattern list files.
    for (ii=0; ii<7; ii++) {
      char patternlist_filename[MAXFILENAME];
      sprintf(patternlist_filename, patternlist_filename_template, ii);
      patf[ii]=openNewPatternFile(patternlist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }
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

      for (ii=0; ii<7; ii++) {
	// Photon list file.
	if (NULL!=plf[ii]) {
	  fits_update_key(plf[ii]->fptr, TDOUBLE, "RA_PNT", &ra,
			  "RA of pointing direction [deg]", &status);
	  fits_update_key(plf[ii]->fptr, TDOUBLE, "DEC_PNT", &dec,
			  "Dec of pointing direction [deg]", &status);
	  fits_update_key(plf[ii]->fptr, TFLOAT, "PA_PNT", &rollangle,
			  "Roll angle [deg]", &status);
	  CHECK_STATUS_BREAK(status);
	}
	
	// Impact list file.
	if (NULL!=ilf[ii]) {
	  fits_update_key(ilf[ii]->fptr, TDOUBLE, "RA_PNT", &ra,
			  "RA of pointing direction [deg]", &status);
	  fits_update_key(ilf[ii]->fptr, TDOUBLE, "DEC_PNT", &dec,
			  "Dec of pointing direction [deg]", &status);
	  fits_update_key(ilf[ii]->fptr, TFLOAT, "PA_PNT", &rollangle,
			  "Roll angle [deg]", &status);
	  CHECK_STATUS_BREAK(status);
	}

	// Event list file.
	fits_update_key(elf[ii]->fptr, TDOUBLE, "RA_PNT", &ra,
			"RA of pointing direction [deg]", &status);
	fits_update_key(elf[ii]->fptr, TDOUBLE, "DEC_PNT", &dec,
			"Dec of pointing direction [deg]", &status);
	fits_update_key(elf[ii]->fptr, TFLOAT, "PA_PNT", &rollangle,
			"Roll angle [deg]", &status);
	CHECK_STATUS_BREAK(status);
	
	// Pattern list file.
	fits_update_key(patf[ii]->fptr, TDOUBLE, "RA_PNT", &ra,
			"RA of pointing direction [deg]", &status);
	fits_update_key(patf[ii]->fptr, TDOUBLE, "DEC_PNT", &dec,
			"Dec of pointing direction [deg]", &status);
	fits_update_key(patf[ii]->fptr, TFLOAT, "PA_PNT", &rollangle,
			"Roll angle [deg]", &status);
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
	
    } else {
      // An explicit attitude file is given.
      for (ii=0; ii<7; ii++) {
	if (NULL!=plf[ii]) {
	  fits_update_key(plf[ii]->fptr, TSTRING, "ATTITUDE", par.Attitude,
			  "attitude file", &status);
	}
	if (NULL!=ilf[ii]) {
	  fits_update_key(ilf[ii]->fptr, TSTRING, "ATTITUDE", par.Attitude,
			  "attitude file", &status);
	}
	fits_update_key(elf[ii]->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
	fits_update_key(patf[ii]->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Timing keywords.
    double buffer_tstop=par.TSTART+par.Exposure;
    double buffer_timezero=0.;
    for (ii=0; ii<7; ii++) {
      // Photon list file.
      if (NULL!=plf[ii]) {
	fits_update_key(plf[ii]->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
			"reference MJD", &status);
	fits_update_key(plf[ii]->fptr, TDOUBLE, "TIMEZERO", &buffer_timezero,
			"time offset", &status);
	fits_update_key(plf[ii]->fptr, TDOUBLE, "TSTART", &par.TSTART,
			"start time", &status);
	fits_update_key(plf[ii]->fptr, TDOUBLE, "TSTOP", &buffer_tstop,
			"stop time", &status);
	CHECK_STATUS_BREAK(status);
      }

      // Impact list file.
      if (NULL!=ilf[ii]) {
	fits_update_key(ilf[ii]->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
			"reference MJD", &status);
	fits_update_key(ilf[ii]->fptr, TDOUBLE, "TIMEZERO", &buffer_timezero,
			"time offset", &status);
	fits_update_key(ilf[ii]->fptr, TDOUBLE, "TSTART", &par.TSTART,
			"start time", &status);
	fits_update_key(ilf[ii]->fptr, TDOUBLE, "TSTOP", &buffer_tstop,
			"stop time", &status);
	CHECK_STATUS_BREAK(status);
      }
      
      // Event list file.
      fits_update_key(elf[ii]->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		      "reference MJD", &status);
      fits_update_key(elf[ii]->fptr, TDOUBLE, "TIMEZERO", &buffer_timezero,
		      "time offset", &status);
      fits_update_key(elf[ii]->fptr, TDOUBLE, "TSTART", &par.TSTART,
		      "start time", &status);
      fits_update_key(elf[ii]->fptr, TDOUBLE, "TSTOP", &buffer_tstop,
		      "stop time", &status);
      CHECK_STATUS_BREAK(status);

      // Pattern list file.
      fits_update_key(patf[ii]->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		      "reference MJD", &status);
      fits_update_key(patf[ii]->fptr, TDOUBLE, "TIMEZERO", &buffer_timezero,
		      "time offset", &status);
      fits_update_key(patf[ii]->fptr, TDOUBLE, "TSTART", &par.TSTART,
		      "start time", &status);
      fits_update_key(patf[ii]->fptr, TDOUBLE, "TSTOP", &buffer_tstop,
		      "stop time", &status);
      CHECK_STATUS_BREAK(status);
    }
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

      // Set the start time for the detector models.
      for (ii=0; ii<7; ii++) {
	setGenDetStartTime(subinst[ii]->det, t0);
      }

      // Loop over photon generation and processing
      // till the time of the photon exceeds the requested
      // time interval.
      do {

	// Photon generation.
	Photon ph;
	int isph=phgen(ac, srccat, MAX_N_SIMPUT, 
		       t0, t1, par.MJDREF, par.dt, 
		       fov7, &ph, &status);
	CHECK_STATUS_BREAK(status);

	// If no photon has been generated, break the loop.
	if (0==isph) break;
	
	// Check if the photon still is within the requested
	// exposre time.
	assert(ph.time<=t1);

	// Randomly assign the photon to one of the 7 sub-telescopes.
	ii=(unsigned int)(sixt_get_random_number(&status)*7);
	CHECK_STATUS_BREAK(status);
	assert(ii<7);

	// If requested, write the photon to the output file.
	if (NULL!=plf[ii]) {
	  status=addPhoton2File(plf[ii], &ph);
	  CHECK_STATUS_BREAK(status);
	}

	// Photon imaging.
	Impact imp;
	int isimg=phimg(subinst[ii]->tel, ac, &ph, &imp, &status);
	CHECK_STATUS_BREAK(status);

	// If the photon is not imaged but lost in the optical system,
	// continue with the next one.
	if (0==isimg) continue;

	// If requested, write the impact to the output file.
	if (NULL!=ilf[ii]) {
	  addImpact2File(ilf[ii], &imp, &status);
	  CHECK_STATUS_BREAK(status);
	}

	// Photon Detection.
	phdetGenDet(subinst[ii]->det, &imp, t1, &status);
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

      // Clear the detectors.
      for (ii=0; ii<7; ii++) {
	phdetGenDet(subinst[ii]->det, NULL, t1, &status);
	CHECK_STATUS_BREAK(status);
	long jj;
	for (jj=0; jj<subinst[ii]->det->pixgrid->ywidth; jj++) {
	  GenDetClearLine(subinst[ii]->det, jj);
	}
      }
      CHECK_STATUS_BREAK(status);

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


    // Use parallel computation via OpenMP.
#pragma omp parallel for reduction(+:status)
    for (ii=0; ii<7; ii++) {
      status=EXIT_SUCCESS;
      // Perform a pattern analysis, only if split events are simulated.
      if (GS_NONE!=subinst[ii]->det->split->type) {
	// Pattern analysis.
	headas_chat(3, "start event pattern analysis ...\n");
	phpat(subinst[ii]->det, elf[ii], patf[ii], par.SkipInvalids, &status);
	//CHECK_STATUS_BREAK(status);
      } else {
	// If no split events are simulated, simply copy the event lists
	// to pattern lists.
	headas_chat(3, "copy events to pattern files ...\n");
	copyEvents2PatternFile(elf[ii], patf[ii], &status);
	//CHECK_STATUS_BREAK(status);
      }
    }
    CHECK_STATUS_BREAK(status);
    
    // Close files in order to save memory.
    for (ii=0; ii<7; ii++) {
      freePhotonListFile(&plf[ii], &status);
      freeImpactListFile(&ilf[ii], &status);
    }

    // Run the event projection.
    headas_chat(3, "start sky projection ...\n");
    for (ii=0; ii<7; ii++) {
      phproj(subinst[ii], ac, patf[ii], par.TSTART, par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  for (ii=0; ii<7; ii++) {
    destroyGenInst    (&subinst[ii], &status);
    destroyPatternFile(&patf[ii],    &status);
    freeEventListFile (&elf[ii],     &status);
    freeImpactListFile(&ilf[ii],     &status);
    freePhotonListFile(&plf[ii],     &status);
  }
  for (ii=0; ii<MAX_N_SIMPUT; ii++) {
    freeSourceCatalog(&srccat[ii], &status);
  }
  freeGTI(&gti);
  freeAttitudeCatalog(&ac);

  if (NULL!=progressfile) {
    fclose(progressfile);
    progressfile=NULL;
  }

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int erosim_getpar(struct Parameters* const par)
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

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile1", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 1");
    return(status);
  } 
  strcpy(par->XMLFile1, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile2", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 2");
    return(status);
  } 
  strcpy(par->XMLFile2, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile3", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 3");
    return(status);
  } 
  strcpy(par->XMLFile3, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile4", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 4");
    return(status);
  } 
  strcpy(par->XMLFile4, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile5", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 5");
    return(status);
  } 
  strcpy(par->XMLFile5, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile6", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 6");
    return(status);
  } 
  strcpy(par->XMLFile6, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile7", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 7");
    return(status);
  } 
  strcpy(par->XMLFile7, sbuffer);
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


