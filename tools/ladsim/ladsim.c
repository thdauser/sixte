#include "ladsim.h"


static void ladphimg(const LAD* const lad,
		     AttitudeCatalog* const ac,
		     PhotonListFile* const plf,
		     LADImpactListFile* const ilf,
		     const double t0,
		     const double exposure,
		     int* const status)
{
  // Calculate the minimum cos-value for sources inside the FOV: 
  // (angle(x0,source) <= 1/2 * diameter)
  const double fov_min_align = cos(lad->fov_diameter/2.); 

  // Scan the entire photon list.  
  while (plf->row < plf->nrows) {

    Photon photon={.time=0.};
      
    // Read an entry from the photon list:
    *status=PhotonListFile_getNextRow(plf, &photon);
    CHECK_STATUS_VOID(*status);

    // Check whether we are still within the requested time interval.
    if (photon.time < t0) continue;
    if (photon.time > t0+exposure) break;

    // Determine telescope pointing direction at the current time.
    struct Telescope telescope;
    telescope.nz = getTelescopeNz(ac, photon.time, status);
    CHECK_STATUS_VOID(*status);

    // Compare the photon direction to the direction of the telescope
    // axis to check whether the photon is inside the FOV.
    Vector photon_direction = unit_vector(photon.ra, photon.dec);
    if (check_fov(&photon_direction, &telescope.nz, fov_min_align)==0) {
      // Photon is inside the FOV!
	
      // Determine telescope data like pointing direction (attitude) etc.
      // The telescope coordinate system consists of an x-, y-, and z-axis.
      getTelescopeAxes(ac, &telescope.nx, &telescope.ny, &telescope.nz, 
		       photon.time, status);
      CHECK_STATUS_BREAK(*status);

      // Apply the geometric vignetting corresponding to the projected
      // detector surface.
      double p = sixt_get_random_number();
      if (p > scalar_product(&photon_direction, &telescope.nz)) {
	// Photon is not detected, since it misses the detector.
	continue;
      }

      // New impact.
      LADImpact impact;
      impact.time = photon.time;
      impact.energy = photon.energy;
      impact.ph_id  = photon.ph_id;
      impact.src_id = photon.src_id;

      // Determine the photon impact position on the detector:
      // Randomly select a panel, module, and element
      impact.panel   = 
	(long)(sixt_get_random_number()*
	       lad->npanels);
      impact.module  = 
	(long)(sixt_get_random_number()*
	       lad->panel[impact.panel]->nmodules);
      impact.element = 	
	(long)(sixt_get_random_number()*
	       lad->panel[impact.panel]->module[impact.module]->nelements);

      // Pointer to the element.
      LADElement* element = 
	lad->panel[impact.panel]->module[impact.module]->element[impact.element];

      // Determine the sensitive area of the element.
      float xwidth = element->xdim - 2.*element->xborder;
      float ywidth = element->ydim - 2.*element->yborder;
      
      // Determine the entrance position into a collimator hole ([m]).
      long col1, row1;
      struct Point2d entrance_position;
      do {
	// Get a random position on the sensitive area of the element.
	entrance_position.x=sixt_get_random_number()*xwidth;
	entrance_position.y=sixt_get_random_number()*ywidth;
	
	// Check if the random position lies within a hole opening.
	LADCollimatorHoleIdx(entrance_position, &col1, &row1);

      } while ((col1<0)||(row1<0));

      // Determine the position on the detector according to the off-axis
      // angle and the orientation of the element ([m]).

      // Determine the photon direction with respect to the telescope
      // coordinate system.
      Vector deviation;
      deviation.x = scalar_product(&telescope.nx, &photon_direction);
      deviation.y = scalar_product(&telescope.ny, &photon_direction);
      deviation.z = scalar_product(&telescope.nz, &photon_direction);

      // Determine the length of the vector to reach from the entrance 
      // position on top of the collimator to the bottom.
      // (The collimator has a thickness of 2 mm.)
      double length = 2.0e-3 / deviation.z;

      // Add the off-axis deviation to the entrance position.
      impact.position.x = 
	entrance_position.x + deviation.x * length;
      impact.position.y = 
	entrance_position.y + deviation.y * length;

      // Check if the impact position still is inside the hole.
      // Otherwise the photon has been absorbed by the walls of the hole.
      // Make sure that it is the SAME hole as before.
      long col2, row2;
      LADCollimatorHoleIdx(impact.position, &col2, &row2);
      if ((col1!=col2)||(row1!=row2)) {
	continue;
      }

      // Check if the impact position still is on the active area of the SDD.
      if ((impact.position.x<0.)||(impact.position.x>xwidth)||
	  (impact.position.y<0.)||(impact.position.y>ywidth)) {
	continue;
      }

      // Insert the impact with the photon data into the list.
      addLADImpact2File(ilf, &impact, status);	    
      CHECK_STATUS_VOID(*status);
      
    } 
    // End of FOV check.
  }
  // END of scanning LOOP over the photon list.
}


static void ladphdet(const LAD* const lad,
		     LADImpactListFile* const ilf,
		     LADSignalListFile* const relf,
		     const double t0,
		     const double exposure,
		     int* const status)
{
  // Loop over all impacts in the FITS file.
  while (ilf->row<ilf->nrows) {

    LADImpact impact;
    getNextLADImpactFromFile(ilf, &impact, status);
    CHECK_STATUS_VOID(*status);

    // Check whether we are still within the requested time interval.
    if (impact.time < t0) continue;
    if (impact.time > t0+exposure) break;


    // Photon detection process according to Campana (2011).
    LADElement* element = 
      lad->panel[impact.panel]->module[impact.module]->element[impact.element];

    // Determine the anode pitch [m].
    float anode_pitch = 
      (element->ydim - 2.*element->yborder)/(element->nanodes-1);

    // Determine the parameters of the charge cloud.
    // Boltzmann constant.
    const double kB = 8.6173324e-5; // [eV/K]
    // Diffusion coefficient.
    //double D = kB * lad->temperature * lad->mobility;
    // Drift velocity.
    double vD = lad->mobility * lad->efield;
    // Charge cloud dispersion.
    const double sigma0 = 20.e-6; // (according to Campana, 2011) [m]
    double sigma = 
      sqrt(2.*kB*lad->temperature*impact.position.x/lad->efield + 
	   pow(sigma0,2.));

    // Drift time.
    double drifttime = impact.position.x / vD;

    // Determine the index of the closest anode strip.
    double y0=impact.position.y/anode_pitch;
    long center_anode=((long)(y0+1)) -1;

    // Determine the measured detector channel (PHA channel) according 
    // to the RMF.
    // The channel is obtained from the RMF using the corresponding
    // HEAdas routine which is based on drawing a random number.
    long channel;
    ReturnChannel(lad->rmf, impact.energy, 1, &channel);

    // Check if the photon is really measured. If the
    // PHA channel returned by the HEAdas RMF function is '-1', 
    // the photon is not detected.
    // This can happen, if the RMF actually is an RSP, i.e. it 
    // includes ARF contributions, e.g., 
    // the detector quantum efficiency and filter transmission.
    if (channel<0) {
      continue; // Photon is not detected.
    }

    // Determine the signal corresponding to the channel according 
    // to the EBOUNDS table.
    float signal = getEBOUNDSEnergy(channel, lad->rmf, 0);
    assert(signal>=0.);

    // Determine which half of the anodes (bottom or top) is affected.
    long min_anode, max_anode;
    if (center_anode < element->nanodes/2) {
      min_anode=0;
      max_anode=element->nanodes/2-1;
    } else {
      min_anode=element->nanodes/2;
      max_anode=element->nanodes-1;
    }

    // Loop over adjacent anodes.
    long ii;
    for (ii =MAX(min_anode,center_anode-2); 
	 ii<=MIN(max_anode,center_anode+2); ii++) {

      LADSignal rev;
      rev.panel   = impact.panel;
      rev.module  = impact.module;
      rev.element = impact.element;
      rev.anode   = ii;
      rev.ph_id[0] = impact.ph_id;
      rev.src_id[0] = impact.src_id;

      // Measured time.
      rev.time = impact.time + drifttime;

      // Measured signal.
      double yi = ii*1.0;
      rev.signal = signal*
	(gaussint(((yi-y0)-0.5)*anode_pitch/sigma)-
	 gaussint(((yi-y0)+0.5)*anode_pitch/sigma));

      // Apply thresholds.
      if (NULL!=lad->threshold_readout_lo_keV) {
	if (rev.signal < *(lad->threshold_readout_lo_keV)) continue;
      }
      if (NULL!=lad->threshold_readout_up_keV) {
	if (rev.signal > *(lad->threshold_readout_up_keV)) continue;
      }

      // Append the new event to the file.
      addLADSignal2File(relf, &rev, status);
      CHECK_STATUS_VOID(*status);
    }
    // END of loop over adjacent anodes.

    // Display program progress.
    if (0==ilf->row%1000) {
      headas_chat(2, "\r %ld of %ld impacts (%.2lf%%) ", 
		  ilf->row, ilf->nrows, ilf->row*100./ilf->nrows);
      fflush(NULL);
    }

  };
  CHECK_STATUS_VOID(*status);
  // END of loop over all impacts in the FITS file.
}


static void ladevents(const LAD* const lad,
		      LADSignalListFile* const relf,
		      LADEventListFile* const elf,
		      int* const status)
{
  // Maximum trigger time between subsequent raw events assigned 
  // to the same event.
  const double dt=1.e-12;

  // List of contributing raw events.
  LADSignal** list=NULL;
  long nlist=0;
  const long maxnlist=10;

  // Error handling loop.
  do { 
    
    // Allocate memory.
    list=(LADSignal**)malloc(maxnlist*sizeof(LADSignal*));
    CHECK_NULL_BREAK(list, *status, "memory allocation for list failed");

    // Loop over all rows in the input file.
    long mm;
    for (mm=0; mm<=relf->nrows; mm++) {
    
      // Read the next raw event from the list.
      LADSignal* rev=NULL;
      if (mm<relf->nrows) {
	rev=getLADSignal(status);
	CHECK_STATUS_BREAK(*status);
	getLADSignalFromFile(relf, mm+1, rev, status);
	CHECK_STATUS_BREAK(*status);
      }

      // Check if the new raw event seems to belong to the same photon event.
      int different=0;
      if ((nlist>0) && (NULL!=rev)) {
	if ((fabs(rev->time-list[0]->time)>dt) ||
	    (rev->panel!=list[0]->panel) ||
	    (rev->module!=list[0]->module) ||
	    (rev->element!=list[0]->element) ||
	    (abs(rev->anode-list[nlist-1]->anode)>1)) {
	  different=1;
	}

	// Check for the different (bottom and top) anode lines.
	long nanodes = 
	  lad->panel[rev->panel]->module[rev->module]->element[rev->element]->nanodes;
	if (((rev->anode>=nanodes/2)&&(list[0]->anode< nanodes/2)) ||
	    ((rev->anode< nanodes/2)&&(list[0]->anode>=nanodes/2))) {
	  different=1;
	}
      }

      // If the new raw event belongs to a different photon event
      // or the end of the input list has been reached, perform 
      // a pattern analysis.
      if ((1==different) || ((mm==relf->nrows)&&(nlist>0))) { 

	// Construct a combined event for output.
	LADEvent* ev = getLADEvent(status);
	CHECK_STATUS_BREAK(*status);
	ev->panel   = list[0]->panel;
	ev->module  = list[0]->module;
	ev->element = list[0]->element;
	ev->time    = list[0]->time;

	long ii;
	long maxidx=0;
	for (ii=0; ii<nlist; ii++) {
	  // Search the anode with the maximum signal.
	  if (list[ii]->signal>list[maxidx]->signal) {
	    maxidx=ii;
	  }
	
	  // Sum the signal contributions.
	  ev->signal += list[ii]->signal;

	  // Set PH_IDs and SRC_IDs.
	  long jj;
	  for (jj=0; jj<NLADSIGNALPHOTONS; jj++) {
	    long kk;
	    for (kk=0; kk<NLADEVENTPHOTONS; kk++) {
	      if (list[ii]->ph_id[jj]==ev->ph_id[kk]) break;
	      if (0==ev->ph_id[kk]) {
		ev->ph_id[kk] =list[ii]->ph_id[jj];
		ev->src_id[kk]=list[ii]->src_id[jj];		
		break;
	      }
	    }
	  }
	}
	ev->anode = list[maxidx]->anode;

	// Add the event to the output file.
	addLADEvent2File(elf, ev, status);
	CHECK_STATUS_BREAK(*status);

	// Release the buffer.
	freeLADEvent(&ev);
	for (ii=0; ii<nlist; ii++) {
	  freeLADSignal(&list[ii]);
	}
	nlist=0;
      }

      // Add the new raw event to the list.
      if (NULL!=rev) {
	if (nlist==maxnlist) {
	  SIXT_ERROR("too many raw events for list buffer");
	  *status=EXIT_FAILURE;
	  break;
	}
	list[nlist]=rev;
	nlist++;
      }
    }
    CHECK_STATUS_VOID(*status);
    // END of loop over all raw events.

  } while(0); // END of error handling loop.
    
  // Release memory.
  if (NULL!=list) {
    long ii;
    for (ii=0; ii<nlist; ii++) {
      if (NULL!=list[ii]) {
	freeLADSignal(&list[ii]);
      }
    }
    free(list);
  }
}


int ladsim_main() 
{
  // Program parameters.
  struct Parameters par;
  
  // Detector setup.
  LAD* lad=NULL;

  // Attitude.
  AttitudeCatalog* ac=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat=NULL;

  // Photon list file.
  PhotonListFile* plf=NULL;

  // Impact list file.
  LADImpactListFile* ilf=NULL;

  // Raw event list file.
  LADSignalListFile* relf=NULL;

  // Recombined event list file.
  LADEventListFile* elf=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("ladsim");
  set_toolversion("0.03");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=ladsim_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Start time for the simulation.
    double t0 = par.TIMEZERO;

    // Determine the appropriate detector XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_LADXMLFile(xml_filename, par.XMLFile);
    CHECK_STATUS_BREAK(status);

    // Determine the photon list output file.
    char ucase_buffer[MAXFILENAME];
    char photonlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.PhotonList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(photonlist_filename, par.OutputStem);
      strcat(photonlist_filename, "_photons.fits");
    } else {
      strcpy(photonlist_filename, par.PhotonList);
    }

    // Determine the impact list output file.
    char impactlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.ImpactList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(impactlist_filename, par.OutputStem);
      strcat(impactlist_filename, "_impacts.fits");
    } else {
      strcpy(impactlist_filename, par.ImpactList);
    }
    
    // Determine the raw event list output file.
    char signallist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.SignalList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(signallist_filename, par.OutputStem);
      strcat(signallist_filename, "_signals.fits");
    } else {
      strcpy(signallist_filename, par.SignalList);
    }

    // Determine the recombined event list output file.
    char eventlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.EventList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(eventlist_filename, par.OutputStem);
      strcat(eventlist_filename, "_events.fits");
    } else {
      strcpy(eventlist_filename, par.EventList);
    }

    // Determine the random number generator seed.
    int seed;
    if (-1!=par.Seed) {
      seed = par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed = (int)time(NULL);
    }

    // Initialize HEADAS random number generator.
    HDmtInit(seed);

    // Load the detector configuration.
    lad=getLADfromXML(xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set up the Attitude.
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a simple pointing attitude.

      // First allocate memory.
      ac=getAttitudeCatalog(&status);
      CHECK_STATUS_BREAK(status);

      ac->entry=(AttitudeEntry*)malloc(2*sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status = EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeCatalog failed");
	break;
      }

      // Set the values of the entries.
      ac->nentries=2;
      ac->entry[0] = defaultAttitudeEntry();
      ac->entry[1] = defaultAttitudeEntry();
      
      ac->entry[0].time = t0;
      ac->entry[1].time = t0+par.Exposure;

      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);
      ac->entry[1].nz = ac->entry[0].nz;

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);
      ac->entry[1].nx = ac->entry[0].nx;

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);
      
      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      if ((ac->entry[0].time > t0) || 
	  (ac->entry[ac->nentries-1].time < t0+par.Exposure)) {
	status=EXIT_FAILURE;
	char msg[MAXMSG];
	sprintf(msg, "attitude data does not cover the "
		"specified period from %lf to %lf!", t0, t0+par.Exposure);
	HD_ERROR_THROW(msg, status);
	break;
      }
    }
    // END of setting up the attitude.

    // Load the SIMPUT X-ray source catalog.
    srccat = loadSourceCatalog(par.Simput, lad->arf, &status);
    CHECK_STATUS_BREAK(status);


#ifdef LAD_OAR
    // Determine the Open Area Ratio of the Collimator on the LAD.
    long kk;
    long try=0, pass=0;
    struct Point2d position;
    float xwidth = 
      lad->panel[0]->module[0]->element[0]->xdim - 
      lad->panel[0]->module[0]->element[0]->xborder*2.;
    float ywidth = 
      lad->panel[0]->module[0]->element[0]->ydim - 
      lad->panel[0]->module[0]->element[0]->yborder*2.;
    for (kk=0; kk<1000000; kk++) {
      long col, row;
      do {
	try++;
	// Get a random position on the element.
	position.x=sixt_get_random_number()*xwidth;
	position.y=sixt_get_random_number()*ywidth;
	
	// Determine the indices of the respective hole.
	LADCollimatorHoleIdx(position, &col, &row);

      } while ((col<0)||(row<0));
      pass++;
    }
    printf("### LAD Open Area Ratio: %lf ###\n", pass*1./try);
#endif


    // --- End of Initialization ---


    // --- Simulation Process ---

    // Open the output photon list file.
    plf=openNewPhotonListFile(photonlist_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);
    fits_update_key(plf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    double dbuffer=0.;
    fits_update_key(plf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		    "time offset", &status);
    CHECK_STATUS_BREAK(status);


    // Photon Generation.
    headas_chat(3, "start photon generation ...\n");
    phgen(ac, srccat, plf, t0, par.Exposure, par.MJDREF, 
	  lad->fov_diameter, &status);
    CHECK_STATUS_BREAK(status);

    // Free the source catalog in order to save memory.
    freeSourceCatalog(&srccat, &status);
    CHECK_STATUS_BREAK(status);

    // Reset internal line counter of photon list file.
    plf->row=0;


    // Open the output impact list file.
    ilf=openNewLADImpactListFile(impactlist_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);
    fits_update_key(ilf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    dbuffer=0.;
    fits_update_key(ilf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		    "time offset", &status);
    CHECK_STATUS_BREAK(status);


    // Photon Imaging.
    headas_chat(3, "start photon imaging ...\n");
    ladphimg(lad, ac, plf, ilf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

    // Close the photon list file in order to save memory.
    freePhotonListFile(&plf, &status);
 

    // Reset internal line counter of impact list file.
    ilf->row=0;


    // Open the output raw event list file.
    relf=openNewLADSignalListFile(signallist_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(relf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);
    fits_update_key(relf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    dbuffer=0.;
    fits_update_key(relf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		    "time offset", &status);
    fits_update_key(relf->fptr, TDOUBLE, "EXPOSURE", &par.Exposure,
		    "exposure time", &status);
    CHECK_STATUS_BREAK(status);

    // Photon Detection.
    headas_chat(3, "start photon detection ...\n");
    ladphdet(lad, ilf, relf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);


    // Close the impact list file in order to save memory.
    freeLADImpactListFile(&ilf, &status);


    // Reset internal line counter of raw event list file.
    relf->row=0;


    // Open the output event list file for recombined events.
    elf=openNewLADEventListFile(eventlist_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(elf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);
    fits_update_key(elf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    dbuffer=0.;
    fits_update_key(elf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		    "time offset", &status);
    fits_update_key(elf->fptr, TDOUBLE, "EXPOSURE", &par.Exposure,
		    "exposure time", &status);
    CHECK_STATUS_BREAK(status);


    // Run the event recombination.
    headas_chat(5, "start event recombination ...\n");
    ladevents(lad, relf, elf, &status);
    CHECK_STATUS_BREAK(status);


    // Run the event projection.
    //headas_chat(5, "start sky projection ...\n");
    // TODO
    //phproj(det, ac, relf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);


    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeLADEventListFile(&elf, &status);
  freeLADSignalListFile(&relf, &status);
  freeLADImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);
  freeSourceCatalog(&srccat, &status);
  freeAttitudeCatalog(&ac);
  freeLAD(&lad);

  // Release HEADAS random number generator:
  HDmtFree();

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}


int ladsim_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("OutputStem", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the filename stem for the output files!\n", 
		   status);
    return(status);
  }
  strcpy(par->OutputStem, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the photon list!\n", status);
    return(status);
  } 
  strcpy(par->PhotonList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("ImpactList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the impact list!\n", status);
    return(status);
  } 
  strcpy(par->ImpactList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("SignalList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the raw event list!\n", status);
    return(status);
  } 
  strcpy(par->SignalList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the event list!\n", status);
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the XML file!\n", status);
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the attitude!\n", status);
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the right ascension of the telescope "
		   "pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the declination of the telescope "
		   "pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the SIMPUT file!\n", status);
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading MJDREF!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("TIMEZERO", &par->TIMEZERO);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading TIMEZERO!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  } 

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the seed for the random number generator!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


