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

#include "athenawfisim.h"


int athenapwfisim_main() 
{

  // Program parameters.
  struct Parameters par;
  
  // Number of detector chips
  unsigned int nchips=4;
  
  unsigned long ii;
  
  // Individual sub-instruments.
  GenInst* subinst[nchips];
  for(ii=0; ii<nchips; ii++){
    subinst[ii]=NULL;
  }

  // Attitude.
  Attitude* ac=NULL;

  // GTI collection.
  GTI* gti=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat[MAX_N_SIMPUT];
  for (ii=0; ii<MAX_N_SIMPUT; ii++) {
    srccat[ii]=NULL;
  }

  // Photon list file.
  PhotonFile* plf=NULL;

  // Impact list file.
  ImpactFile* ilf=NULL;

  // Event list files.
  EventFile* elf[nchips];
  for(ii=0; ii<nchips; ii++){
    elf[ii]=NULL;
  }

  // Pattern list files.
  EventFile* patf[nchips];
  for(ii=0; ii<nchips; ii++){
    patf[ii]=NULL;
  }

  // Output file for progress status.
  FILE* progressfile=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("athenapwfisim");
  set_toolversion("0.08");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=athenapwfisim_getpar(&par);
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
    char rawdata_filename_template[MAXFILENAME];
    strcpy(ucase_buffer, par.RawData);
    strtoupper(ucase_buffer);
    strcpy(rawdata_filename_template, par.Prefix);
    strcat(rawdata_filename_template, "chip%d_");
    int delete_rawdata = 0;
    if (0==strcmp(ucase_buffer,"NONE")) {
    	delete_rawdata = 1;
      strcat(rawdata_filename_template, "raw.fits");
    } else {
      strcat(rawdata_filename_template, par.RawData);
    }

    // Determine the pattern list output file.
    char evtfile_filename_template[MAXFILENAME];
    strcpy(ucase_buffer, par.EvtFile);
    strtoupper(ucase_buffer);
    strcpy(evtfile_filename_template, par.Prefix);
    strcat(evtfile_filename_template, "chip%d_");
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcat(evtfile_filename_template, "evt.fits");
    } else {
      strcat(evtfile_filename_template, par.EvtFile);
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

    // Load the configurations of all five chips from the XML 
    // definition files.
    for (ii=0; ii<nchips; ii++) {
      // Check if a particular XML file is given for this 
      // chip. If not, use the default XML file.
      char buffer[MAXFILENAME];
      char ubuffer[MAXFILENAME];
      char default_filename[MAXFILENAME];
      switch (ii) {
      case 0:
	strcpy(buffer, par.XMLFile0);
	strcpy(default_filename, "depfet_b_1l_ff_chip0.xml");
	break;
      case 1:
	strcpy(buffer, par.XMLFile1);
	strcpy(default_filename, "depfet_b_1l_ff_chip1.xml");
	break;
      case 2:
	strcpy(buffer, par.XMLFile2);
	strcpy(default_filename, "depfet_b_1l_ff_chip2.xml");
	break;
      case 3:
	strcpy(buffer, par.XMLFile3);
	strcpy(default_filename, "depfet_b_1l_ff_chip3.xml");
	break;
      default:
	break;
      }
      strcpy(ubuffer, buffer);
      strtoupper(ubuffer);
      if (0==strcmp(ubuffer, "NONE")) {
	strcpy(buffer, SIXT_DATA_PATH);
	strcat(buffer, "/instruments/athena/wfi_wo_filter/");
	strcat(buffer, default_filename);
      }

      // Load the instrument configuration either with the
      // specific (if available) or the default XML file.
      subinst[ii]=loadGenInst(buffer, seed, &status);
      CHECK_STATUS_BREAK(status);

      // Set the usage of the detector background according to
      // the respective program parameter.
      setGenDetIgnoreBkg(subinst[ii]->det, !par.Background);
    }
    CHECK_STATUS_BREAK(status);
    

    // Set up the Attitude.
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if ((strlen(par.Attitude)==0)||(0==strcmp(ucase_buffer, "NONE"))) {
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

    // Get a GTI.
    gti=getGTIFromFileOrContinuous(par.GTIfile, 
				   par.TSTART, par.TSTART+par.Exposure,
				   par.MJDREF, &status);
    CHECK_STATUS_BREAK(status);

    // Load the SIMPUT X-ray source catalogs.
    srccat[0]=loadSourceCatalog(par.Simput, subinst[0]->tel->arf, &status);
    CHECK_STATUS_BREAK(status);

    // Optional 2nd catalog.
    if (strlen(par.Simput2)>0) {
      strcpy(ucase_buffer, par.Simput2);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[1]=loadSourceCatalog(par.Simput2, subinst[0]->tel->arf, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 3rd catalog.
    if (strlen(par.Simput3)>0) {
      strcpy(ucase_buffer, par.Simput3);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
	srccat[2]=loadSourceCatalog(par.Simput3, subinst[0]->tel->arf, &status);
	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 4th catalog.
    if (strlen(par.Simput4)>0) {
      strcpy(ucase_buffer, par.Simput4);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[3]=loadSourceCatalog(par.Simput4, subinst[0]->tel->arf, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 5th catalog.
    if (strlen(par.Simput5)>0) {
      strcpy(ucase_buffer, par.Simput5);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[4]=loadSourceCatalog(par.Simput5, subinst[0]->tel->arf, &status);
      	CHECK_STATUS_BREAK(status);
      }
    }

    // Optional 6th catalog.
    if (strlen(par.Simput6)>0) {
      strcpy(ucase_buffer, par.Simput6);
      strtoupper(ucase_buffer);
      if (0!=strcmp(ucase_buffer, "NONE")) {
      	srccat[5]=loadSourceCatalog(par.Simput6, subinst[0]->tel->arf, &status);
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

    double tstop=gti->stop[gti->ngti-1];

    // Open the output photon list files.
    if (strlen(photonlist_filename)>0) {
      plf=openNewPhotonFile(photonlist_filename,
			    telescop, instrume, "Normal", 
			    subinst[0]->tel->arf_filename,
			    subinst[0]->det->rmf_filename,
			    par.MJDREF, 0.0, par.TSTART, tstop,
			    par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list files.
    if (strlen(impactlist_filename)>0) {
      ilf=openNewImpactFile(impactlist_filename, 
			    telescop, instrume, "Normal", 
			    subinst[0]->tel->arf_filename,
			    subinst[0]->det->rmf_filename,
			    par.MJDREF, 0.0, par.TSTART, tstop,
			    par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list files.
    for (ii=0; ii<nchips; ii++) {
      char rawdata_filename[MAXFILENAME];
      sprintf(rawdata_filename, rawdata_filename_template, ii);
      elf[ii]=openNewEventFile(rawdata_filename,
			       telescop, instrume, "Normal", 
			       subinst[0]->tel->arf_filename,
			       subinst[0]->det->rmf_filename,
			       par.MJDREF, 0.0, par.TSTART, tstop,
			       subinst[ii]->det->pixgrid->xwidth,
			       subinst[ii]->det->pixgrid->ywidth,
			       par.clobber, &status);
      CHECK_STATUS_BREAK(status);

      // Define the event list file as output file for the respective
      // detector chip.
      setGenDetEventFile(subinst[ii]->det, elf[ii]);
    }
    CHECK_STATUS_BREAK(status);

    // Open the output pattern list files.
    for (ii=0; ii<nchips; ii++) {
      char evtfile_filename[MAXFILENAME];
      sprintf(evtfile_filename, evtfile_filename_template, ii);
      patf[ii]=openNewEventFile(evtfile_filename,
				telescop, instrume, "Normal", 
				subinst[0]->tel->arf_filename,
				subinst[0]->det->rmf_filename,
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

      for (ii=0; ii<nchips; ii++) {
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

      if (NULL!=plf) {
	fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
      }
      if (NULL!=ilf) {
	fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
      }

      for (ii=0; ii<nchips; ii++) {
	fits_update_key(elf[ii]->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
	fits_update_key(patf[ii]->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
    }

    // Event type.
    for (ii=0; ii<nchips; ii++) {
      fits_update_key(elf[ii]->fptr, TSTRING, "EVTYPE", "PIXEL",
		      "event type", &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);  

    // TLMIN and TLMAX of PI column.
    for (ii=0; ii<nchips; ii++) {
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
    for (ii=0; ii<nchips; ii++) {
      fits_update_key(patf[ii]->fptr, TDOUBLE, "EXPOSURE", &totalsimtime,
		      "exposure time [s]", &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Loop over all intervals in the GTI collection.
    double simtime=0.;
    int gtibin=0;
    do {
      // Currently regarded interval.
      double t0=gti->start[gtibin];
      double t1=gti->stop[gtibin];

      // Set the start time for the detector models.
      for (ii=0; ii<nchips; ii++) {
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
		       subinst[0]->tel->fov_diameter, &ph, &status);
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
	int isimg=phimg(subinst[0]->tel, ac, &ph, &imp, &status);
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
	for (ii=0; ii<nchips; ii++) {
	  phdetGenDet(subinst[ii]->det, &imp, t1, &status);
	  CHECK_STATUS_BREAK(status);
	}
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

      // Clear the detectors.
      for (ii=0; ii<nchips; ii++) {
	phdetGenDet(subinst[ii]->det, NULL, t1, &status);
	CHECK_STATUS_BREAK(status);
	long jj;
	for (jj=0; jj<subinst[ii]->det->pixgrid->ywidth; jj++) {
	  GenDetClearLine(subinst[ii]->det, jj);
	}
      }
      CHECK_STATUS_BREAK(status);

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


    // Use parallel computation via OpenMP.
#pragma omp parallel for reduction(+:status)
    for (ii=0; ii<nchips; ii++) {
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
	copyEventFile(elf[ii], patf[ii],
		      subinst[ii]->det->threshold_event_lo_keV,
		      subinst[ii]->det->threshold_pattern_up_keV,
		      &status);
	//CHECK_STATUS_BREAK(status);
	fits_update_key(patf[ii]->fptr, TSTRING, "EVTYPE", "PATTERN", 
			"event type", &status);
	//CHECK_STATUS_BREAK(status);
      }
      //CHECK_STATUS_BREAK(status);
      // END of loop over all events in the list.
    }
    CHECK_STATUS_BREAK(status);
    
    // Store the GTI extension in the event files.
    for (ii=0; ii<nchips; ii++) {
      saveGTIExt(elf[ii]->fptr, "STDGTI", gti, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Close files in order to save memory.
    freePhotonFile(&plf, &status);
    freeImpactFile(&ilf, &status);
    CHECK_STATUS_BREAK(status);
    for (ii=0; ii<nchips; ii++) {
      freeEventFile(&elf[ii], &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Run the event projection.
    headas_chat(3, "start sky projection ...\n");
    for (ii=0; ii<nchips; ii++) {
      phproj(subinst[ii], ac, patf[ii], par.TSTART, par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Store the GTI extension in the pattern files.
    for (ii=0; ii<nchips; ii++) {
      saveGTIExt(patf[ii]->fptr, "STDGTI", gti, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // --- End of simulation process ---

    // remove RawData files if not requested
    if (delete_rawdata){
    	for (ii=0; ii<nchips; ii++) {
    		char rawdata_filename[MAXFILENAME];
    		sprintf(rawdata_filename, rawdata_filename_template, ii);
    		headas_chat(5,"removing unwanted RawData file %s \n",rawdata_filename);
    		status = remove (rawdata_filename);
    		CHECK_STATUS_BREAK(status);
    	}
    }


  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");


  // Release memory.
  freeImpactFile(&ilf, &status);
  freePhotonFile(&plf, &status);
  for (ii=0; ii<nchips; ii++) {
    destroyGenInst(&subinst[ii], &status);
    freeEventFile(&patf[ii], &status);
    freeEventFile(&elf[ii], &status);
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

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int athenapwfisim_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // check if any obsolete keywords are given
  sixt_check_obsolete_keyword(&status);
  CHECK_STATUS_RET(status,EXIT_FAILURE);

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
    SIXT_ERROR("failed reading the name of the raw data output file");
    return(status);
  }
  strcpy(par->RawData, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output event file");
    return(status);
  }
  strcpy(par->EvtFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile0", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file 0");
    return(status);
  }
  strcpy(par->XMLFile0, sbuffer);
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

