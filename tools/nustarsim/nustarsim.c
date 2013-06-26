#include "nustarsim.h"


int nustarsim_main() 
{
  // Program parameters.
  struct Parameters par;
  
  // Individual sub-instruments.
  GenInst* subinst[2]={NULL, NULL};

  // Fake telescope ARF with the composite effective area of both
  // sub-telescopes. This is used for the photon generation
  // from the SIMPUT catalog.
  struct ARF* arf2=NULL;
  
  // FoV of the two telescopes combined. If there is not misalignment,
  // it corresponds to the FoV of a single sub-telescope.
  float fov2;

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

  // Photon list files.
  PhotonListFile* plf[2]={NULL, NULL};

  // Impact list file.
  ImpactListFile* ilf[2]={NULL, NULL};

  // Event list file.
  EventListFile* elf[2]={NULL, NULL};

  // Pattern list file.
  PatternFile* patf[2]={NULL, NULL};

  // Output file for progress status.
  FILE* progressfile=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("nustarsim");
  set_toolversion("0.04");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=nustarsim_getpar(&par);
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
    strcat(eventlist_filename_template, "module%d_");
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
    strcat(patternlist_filename_template, "module%d_");
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
		     "NUSTAR", "NUSTAR", "", &status);
    CHECK_STATUS_BREAK(status);

    // Load the configurations of both sub-instruments.
    for (ii=0; ii<2; ii++) {
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

      // Make sure that no split model is selected.
      // This option is currently not supported for the NuSTAR
      // instrument model.
      assert(GS_NONE==subinst[ii]->det->split->type);
    }
    CHECK_STATUS_BREAK(status);
    
    // Determine a fake ARF with the combined effective area of
    // both sub-telescopes.
    arf2=(struct ARF*)malloc(sizeof(struct ARF));
    CHECK_NULL_BREAK(arf2, status, "memory allocation for fake ARF failed");
    arf2->NumberEnergyBins=subinst[0]->tel->arf->NumberEnergyBins;
    arf2->LowEnergy= (float*)malloc(arf2->NumberEnergyBins*sizeof(float));
    arf2->HighEnergy=(float*)malloc(arf2->NumberEnergyBins*sizeof(float));
    arf2->EffArea=   (float*)malloc(arf2->NumberEnergyBins*sizeof(float));
    long kk;
    for (kk=0; kk<arf2->NumberEnergyBins; kk++) {
      arf2->LowEnergy[kk] =subinst[0]->tel->arf->LowEnergy[kk];
      arf2->HighEnergy[kk]=subinst[0]->tel->arf->HighEnergy[kk];
      arf2->EffArea[kk]   =subinst[0]->tel->arf->EffArea[kk] *2.0; // Factor 2 !
    }
    strcpy(arf2->ARFVersion, subinst[0]->tel->arf->ARFVersion);
    strcpy(arf2->Telescope , subinst[0]->tel->arf->Telescope );
    strcpy(arf2->Instrument, subinst[0]->tel->arf->Instrument);
    strcpy(arf2->Detector  , subinst[0]->tel->arf->Detector  );
    strcpy(arf2->Filter    , subinst[0]->tel->arf->Filter    );
    strcpy(arf2->ARFExtensionName, subinst[0]->tel->arf->ARFExtensionName);
    
    // The FoV is the same as for an individual sub-telescope.
    // We ignore any misalignment for the moment.
    fov2=subinst[0]->tel->fov_diameter;
    
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

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);

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
    srccat[0]=loadSourceCatalog(par.Simput, arf2, &status);
    CHECK_STATUS_BREAK(status);

    // Optional 2nd catalog.
    if (strlen(par.Simput2)>0) {
      strcpy(ucase_buffer, par.Simput2);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[1]=loadSourceCatalog(par.Simput2, arf2, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 3rd catalog.
    if (strlen(par.Simput3)>0) {
      strcpy(ucase_buffer, par.Simput3);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[2]=loadSourceCatalog(par.Simput3, arf2, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 4th catalog.
    if (strlen(par.Simput4)>0) {
      strcpy(ucase_buffer, par.Simput4);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[3]=loadSourceCatalog(par.Simput4, arf2, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 5th catalog.
    if (strlen(par.Simput5)>0) {
      strcpy(ucase_buffer, par.Simput5);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[4]=loadSourceCatalog(par.Simput5, arf2, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 6th catalog.
    if (strlen(par.Simput6)>0) {
      strcpy(ucase_buffer, par.Simput6);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[5]=loadSourceCatalog(par.Simput6, arf2, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // --- End of Initialization ---


    // --- Open and set up files ---

    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=subinst[0]->telescop) {
      strcpy(telescop, subinst[0]->telescop);
    }
    if (NULL!=subinst[0]->instrume) {
      strcpy(instrume, subinst[0]->instrume);
    }
    
    double tstop;
    if (NULL==gti) {
      tstop=par.TSTART+par.Exposure;
    } else {
      tstop=gti->stop[gti->nentries-1];
    }

    // Open the output photon list files.
    if (strlen(photonlist_filename_template)>0) {
      for (ii=0; ii<2; ii++) {
	char photonlist_filename[MAXFILENAME];
	sprintf(photonlist_filename, photonlist_filename_template, ii);
	plf[ii]=openNewPhotonListFile(photonlist_filename, 
				      telescop, instrume, "Normal", 
				      par.MJDREF, 0.0, par.TSTART, tstop,
				      par.clobber, &status);
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list files.
    if (strlen(impactlist_filename_template)>0) {
      for (ii=0; ii<2; ii++) {
	char impactlist_filename[MAXFILENAME];
	sprintf(impactlist_filename, impactlist_filename_template, ii);
	ilf[ii]=openNewImpactListFile(impactlist_filename, 
				      telescop, instrume, "Normal", 
				      par.MJDREF, 0.0, par.TSTART, tstop,
				      par.clobber, &status);
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list files.
    for (ii=0; ii<2; ii++) {
      char eventlist_filename[MAXFILENAME];
      sprintf(eventlist_filename, eventlist_filename_template, ii);
      elf[ii]=openNewEventListFile(eventlist_filename, 
				   telescop, instrume, "Normal", 
				   par.MJDREF, 0.0, par.TSTART, tstop,
				   subinst[ii]->det->pixgrid->xwidth,
				   subinst[ii]->det->pixgrid->ywidth,
				   par.clobber, &status);
      CHECK_STATUS_BREAK(status);

      // Define the event list file as output file for the respective
      // detector.
      setGenDetEventListFile(subinst[ii]->det, elf[ii]);
    }
    CHECK_STATUS_BREAK(status);

    // Open the output pattern list files.
    for (ii=0; ii<2; ii++) {
      char patternlist_filename[MAXFILENAME];
      sprintf(patternlist_filename, patternlist_filename_template, ii);
      patf[ii]=openNewPatternFile(patternlist_filename, 
				  telescop, instrume, "Normal",
				  par.MJDREF, 0.0, par.TSTART, tstop,
				  subinst[ii]->det->pixgrid->xwidth,
				  subinst[ii]->det->pixgrid->ywidth,
				  par.clobber, &status);
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
      ra *=180./M_PI;
      dec*=180./M_PI;
      rollangle*=180./M_PI;

      for (ii=0; ii<2; ii++) {
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
      for (ii=0; ii<2; ii++) {
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

    // TLMIN and TLMAX of PI column.
    for (ii=0; ii<2; ii++) {
      char keystr[MAXMSG];
      long value;
      sprintf(keystr, "TLMIN%d", elf[ii]->cpi);
      value=subinst[ii]->det->rmf->FirstChannel;
      fits_update_key(elf[ii]->fptr, TLONG, keystr, &value, "", &status);
      sprintf(keystr, "TLMAX%d", elf[ii]->cpi);
      value=subinst[ii]->det->rmf->FirstChannel+subinst[ii]->det->rmf->NumberChannels-1;
      fits_update_key(elf[ii]->fptr, TLONG, keystr, &value, "", &status);
      CHECK_STATUS_BREAK(status);
    
      sprintf(keystr, "TLMIN%d", patf[ii]->cpi);
      value=subinst[ii]->det->rmf->FirstChannel;
      fits_update_key(patf[ii]->fptr, TLONG, keystr, &value, "", &status);
      sprintf(keystr, "TLMAX%d", patf[ii]->cpi);
      value=subinst[ii]->det->rmf->FirstChannel+subinst[ii]->det->rmf->NumberChannels-1;
      fits_update_key(patf[ii]->fptr, TLONG, keystr, &value, "", &status);
      CHECK_STATUS_BREAK(status);  
    }

    // Timing keywords.
    double buffer_tstop=par.TSTART+par.Exposure;
    double buffer_timezero=0.;
    for (ii=0; ii<2; ii++) {
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
	t1=par.TSTART+par.Exposure;
      } else {
	t0=gti->start[gtibin];
	t1=gti->stop[gtibin];
      }

      // Set the start time for the detector models.
      for (ii=0; ii<2; ii++) {
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
		       fov2, &ph, &status);
	CHECK_STATUS_BREAK(status);

	// If no photon has been generated, break the loop.
	if (0==isph) break;
	
	// Check if the photon still is within the requested
	// exposre time.
	assert(ph.time<=t1);

	// Randomly assign the photon to one of both sub-telescopes.
	ii=(unsigned int)(sixt_get_random_number(&status)*2.0);
	CHECK_STATUS_BREAK(status);
	assert(ii<2);

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

	// Skip if the impact position is located in the gap of 
	// the 2x2 module. Note that this approach does not work
	// properly if split events between neighboring pixels are
	// enabled. In this case there needs to be a particular
	// check in the detection routine on the affected pixels.
	const double detector_offset=3300.e-6;
	if (((imp.position.x>=-300.e-6-detector_offset)&&
	     (imp.position.x<=300.e-6-detector_offset))||
	    ((imp.position.y>=-300.e-6+detector_offset)&&
	     (imp.position.y<=300.e-6+detector_offset))) {
	  continue;
	}

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
      for (ii=0; ii<2; ii++) {
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
    for (ii=0; ii<2; ii++) {
      status=EXIT_SUCCESS;

      // Apply dead time, charge pump reset, and shield anticoincidence 
      // intervals, while producing a pattern list from the event list. 
      // The event list contains all events neglecting dead time, charge 
      // pump resets, and shield vetos, while the pattern list contains
      // only events after the dead time application.
      headas_chat(3, "apply dead time ...\n");
      double last_time=0.;
      double veto_time=0.;
      const double veto_interval=500.e-6;
      const double veto_rate=28.;

      // Loop over all rows in the event file.
      long row;
      for (row=0; row<elf[ii]->nrows; row++) {
	// Buffers.
	Event event;
	Pattern pattern;
	
	// Read an event from the input list.
	getEventFromFile(elf[ii], row+1, &event, &status);
	CHECK_STATUS_BREAK(status);
    
	// Make sure that the event is not located in the gap
	// between the 2x2 hybrids.
	assert((event.rawx<160)||(event.rawx>164));
	assert((event.rawy<160)||(event.rawy>164));
	
	// Check if the event falls within an interval of charge
	// pump reset. Resets take place every millisecond and 
	// last for 20mus.
	double dt=event.time - ((long)(event.time*1000.))*0.001;
	if (dt<0.02e-3) continue;


	// Apply the shield veto time.
	while (event.time-veto_time>veto_interval) {
	  veto_time+=rndexp(1./veto_rate, &status);
	  CHECK_STATUS_BREAK(status);
	}
	CHECK_STATUS_BREAK(status);
	
	if ((event.time-veto_time>=0.) &&
	    (event.time-veto_time<=veto_interval)) {
	  continue;
	}
	

	// Apply the dead time (event processing time).
	if (0==row) {
	  last_time=event.time;
	} else if (event.time-last_time<=2.5e-3) {
	  continue;
	}

	// The event is detected.
	last_time=event.time;

	// Copy event data to pattern.
	pattern.rawx   =event.rawx;
	pattern.rawy   =event.rawy;
	pattern.time   =event.time;
	pattern.frame  =event.frame;
	pattern.pi     =event.pi;
	pattern.signal =event.signal;
	pattern.ra     =0.;
	pattern.dec    =0.;
	pattern.npixels=1;
	pattern.type   =0;
    
	pattern.pileup =0;
	int jj;
	for (jj=0; (jj<NEVENTPHOTONS)&&(jj<NPATTERNPHOTONS); jj++){
	  pattern.ph_id[jj] =event.ph_id[jj];
	  pattern.src_id[jj]=event.src_id[jj];
	  
	  if ((jj>0)&&(pattern.ph_id[jj]!=0)) {
	    pattern.pileup=1;
	  }
	}
	    
	pattern.signals[0]=0.;
	pattern.signals[1]=0.;
	pattern.signals[2]=0.;
	pattern.signals[3]=0.;
	pattern.signals[4]=event.signal;
	pattern.signals[5]=0.;
	pattern.signals[6]=0.;
	pattern.signals[7]=0.;
	pattern.signals[8]=0.;
	
	// Add the new pattern to the output file.
	addPattern2File(patf[ii], &pattern, &status);	  
	CHECK_STATUS_BREAK(status);
      }
      //CHECK_STATUS_BREAK(status);
      // END of loop over all events in the list.
    }
    CHECK_STATUS_BREAK(status);
    
    // Close files in order to save memory.
    for (ii=0; ii<2; ii++) {
      freePhotonListFile(&plf[ii], &status);
      freeImpactListFile(&ilf[ii], &status);
    }

    // Run the event projection.
    headas_chat(3, "start sky projection ...\n");
    for (ii=0; ii<2; ii++) {
      phproj(subinst[ii], ac, patf[ii], par.TSTART, par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  for (ii=0; ii<2; ii++) {
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
  freeAttitude(&ac);

  if (NULL!=progressfile) {
    fclose(progressfile);
    progressfile=NULL;
  }

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int nustarsim_getpar(struct Parameters* const par)
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

