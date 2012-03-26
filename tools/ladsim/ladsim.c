#include "ladsim.h"


struct LADSignalListItem {
  /** Signal entry. */
  LADSignal* signal;
  
  /** Pointer to the next item. */
  struct LADSignalListItem* next;
};


static inline int ladphimg(const LAD* const lad,
			   AttitudeCatalog* const ac,
			   Photon* const ph,
			   LADImpact* const imp,
			   int* const status)
{
  // Calculate the minimum cos-value for sources inside the FOV: 
  // (angle(x0,source) <= 1/2 * diameter)
  const double fov_min_align = cos(lad->fov_diameter/2.); 

  // Determine telescope pointing direction at the current time.
  struct Telescope telescope;
  telescope.nz = getTelescopeNz(ac, ph->time, status);
  CHECK_STATUS_RET(*status, 0);

  // Compare the photon direction to the direction of the telescope
  // axis to check whether the photon is inside the FOV.
  Vector photon_direction = unit_vector(ph->ra, ph->dec);
  if (check_fov(&photon_direction, &telescope.nz, fov_min_align)==0) {
    // Photon is inside the FOV!
	
    // Determine telescope data like pointing direction (attitude) etc.
    // The telescope coordinate system consists of an x-, y-, and z-axis.
    getTelescopeAxes(ac, &telescope.nx, &telescope.ny, &telescope.nz, 
		     ph->time, status);
    CHECK_STATUS_RET(*status, 0);

    // Apply the geometric vignetting corresponding to the projected
    // detector surface.
    double p = sixt_get_random_number();
    if (p > scalar_product(&photon_direction, &telescope.nz)) {
      // Photon is not detected, since it misses the detector.
      return(0);
    }

    // New impact.
    imp->time = ph->time;
    imp->energy = ph->energy;
    imp->ph_id  = ph->ph_id;
    imp->src_id = ph->src_id;
    
    // Determine the photon impact position on the detector:
    // Randomly select a panel, module, and element
    imp->panel   = 
      (long)(sixt_get_random_number()*lad->npanels);
    imp->module  = 
      (long)(sixt_get_random_number()*lad->panel[imp->panel]->nmodules);
    imp->element = 	
      (long)(sixt_get_random_number()*
	     lad->panel[imp->panel]->module[imp->module]->nelements);
    
    // Pointer to the element.
    LADElement* element = 
      lad->panel[imp->panel]->module[imp->module]->element[imp->element];
    
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
    imp->position.x = 
      entrance_position.x + deviation.x * length;
    imp->position.y = 
      entrance_position.y + deviation.y * length;
    
    // Check if the impact position still is inside the hole.
    // Otherwise the photon has been absorbed by the walls of the hole.
    // Make sure that it is the SAME hole as before.
    long col2, row2;
    LADCollimatorHoleIdx(imp->position, &col2, &row2);
    if ((col1!=col2)||(row1!=row2)) {
      return(0);
    }

    // Check if the impact position still is on the active area of the SDD.
    if ((imp->position.x<0.)||(imp->position.x>xwidth)||
	(imp->position.y<0.)||(imp->position.y>ywidth)) {
      return(0);
    }
  } 
  // End of FOV check.

  return(1);
}


static inline int ladphdet(const LAD* const lad,
			   LADImpact* const imp,
			   LADSignal* const signals)
{
  // Photon detection process according to Campana (2011).
  LADElement* element = 
    lad->panel[imp->panel]->module[imp->module]->element[imp->element];

  // Determine the anode pitch [m].
  float anode_pitch = 
    (element->ydim - 2.*element->yborder)/(element->nanodes/2);

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
    sqrt(2.*kB*lad->temperature*imp->position.x/lad->efield + 
	 pow(sigma0,2.));

  // Drift time.
  double drifttime = imp->position.x / vD;

  // Determine the index of the closest anode strip.
  float xwidth = element->xdim - 2.*element->xborder;
  double y0;
  if (imp->position.x<0.5*xwidth) {
    y0=imp->position.y/anode_pitch;
  } else {
    y0=imp->position.y/anode_pitch + element->nanodes/2;
  }
  long center_anode=((long)(y0+1)) -1;

  // Determine the measured detector channel (PHA channel) according 
  // to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  returnRMFChannel(lad->rmf, imp->energy, &channel);

  // Check if the photon is really measured. If the
  // PHA channel returned by the HEAdas RMF function is '-1', 
  // the photon is not detected.
  // This can happen, if the RMF actually is an RSP, i.e. it 
  // includes ARF contributions, e.g., 
  // the detector quantum efficiency and filter transmission.
  if (channel<0) {
    headas_chat(5, "### undetected photon\n");
    return(0); // Photon is not detected.
  } else if (channel<lad->rmf->FirstChannel) {
    headas_printf("### wrong channel number !\n");
    return(0);
  }

  // Determine the signal corresponding to the channel according 
  // to the EBOUNDS table.
  float signal = getEBOUNDSEnergy(channel, lad->rmf, 0);
  assert(signal>=0.);

  // Determine which half of the anodes (bottom or top) is affected.
  long min_anode, max_anode;
  if (center_anode < element->nanodes/2) {
    min_anode=MAX(0                   , center_anode-2);
    max_anode=MIN(element->nanodes/2-1, center_anode+2);
  } else {
    min_anode=MAX(element->nanodes/2, center_anode-2);
    max_anode=MIN(element->nanodes-1, center_anode+2);
  }
  int n_anodes=max_anode-min_anode+1;
  assert(n_anodes<=5);
    
  // Loop over adjacent anodes.
  int ii, jj=0; // (lies within [0,4])
  for (ii=0; ii<n_anodes; ii++) {

    // Determine the signal fraction at this anode.
    double yi=(ii+min_anode)*1.0;
    double fraction=0.;
    if (ii>0) {
      fraction =gaussint((yi-y0)*anode_pitch/sigma);
    } else {
      fraction =1.;
    }
    if (ii<n_anodes-1) {
      fraction -=gaussint((yi-y0+1.0)*anode_pitch/sigma);
    }
    signals[jj].signal = fraction*signal;

    // Apply thresholds.
    if (NULL!=lad->threshold_readout_lo_keV) {
      if (signals[jj].signal < *(lad->threshold_readout_lo_keV)) {
	continue;
      }
    }
    if (NULL!=lad->threshold_readout_up_keV) {
      if (signals[jj].signal > *(lad->threshold_readout_up_keV)) {
	continue;
      }
    }

    // Determine the point of time at which the charge is 
    // collected on the anode.
    signals[jj].time = imp->time + drifttime;

    signals[jj].panel   = imp->panel;
    signals[jj].module  = imp->module;
    signals[jj].element = imp->element;
    signals[jj].anode   = ii+min_anode;
    signals[jj].ph_id[0]  = imp->ph_id;
    signals[jj].src_id[0] = imp->src_id;
    int kk;
    for (kk=1; kk<NLADSIGNALPHOTONS; kk++) {
      signals[jj].ph_id[kk] = 0;
      signals[jj].src_id[kk]= 0;
    }

    // Increment counter for signal buffer.
    jj++;

  }
  // END of loop over adjacent anodes.

  // If no signal is above the threshold, we can skip here.
  if (jj==0) return(0);


  // Apply the ASIC dead time.
  // Determine the ASIC and the pin of the ASIC to which the 
  // central anode is attached.
  int asic=(int)(center_anode/lad->asic_channels);
  int pin =      center_anode%lad->asic_channels;

  // Check if the bin is at the border of the ASIC and whether
  // the neighboring ASIC has to be read out, too.
  int asic2=-1;
  if ((pin<2) && (asic!=0) && (asic!=element->nasics/2)) {
    asic2 = asic-1;
  } else if ((pin>lad->asic_channels-3) && 
	     (asic!=element->nasics/2-1) &&
	     (asic!=element->nasics-1)) {
    asic2 = asic+1;
  }

  // Check if the event happens after the coincidence time, but
  // during the dead time.
  if (element->asic_readout_time[asic]>0.) {
    if ((signals[0].time-element->asic_readout_time[asic]>lad->coincidencetime) &&
	(signals[0].time-element->asic_readout_time[asic]<
	 lad->coincidencetime+lad->deadtime)) {
      return(0);
    }
  } else if (asic2>=0) {
    if (element->asic_readout_time[asic2]>0.) {
      if ((signals[0].time-element->asic_readout_time[asic2]>lad->coincidencetime) &&
	  (signals[0].time-element->asic_readout_time[asic2]<
	   lad->coincidencetime+lad->deadtime)) {
	return(0);
      }
    }
  }

  // Set the time of this ASIC readout.
  if ((element->asic_readout_time[asic]==0.) ||
      (signals[0].time-element->asic_readout_time[asic]>lad->coincidencetime)) {
    element->asic_readout_time[asic]=signals[0].time;
    if (asic2>=0) {
      element->asic_readout_time[asic2]=signals[0].time;
    }
  }
  // END of dead time application.

  return(jj);
}


/** Return recombined events from the detected signals. */
static inline int ladevrecomb(const LAD* const lad,
			      LADSignal* const signal,
			      LADEvent* const event,
			      int* const status)
{
  // List of raw event signals.
  static struct LADSignalListItem* first=NULL;

  // Flag if an event is complete (there will be no further
  // signal contributions).
  int complete=0;
  if (NULL!=first) {
    if (NULL==signal) {
      complete=1;
    } else if (signal->time-first->signal->time>lad->coincidencetime) {
      complete=1;
    }
  }

  // Produce a new event.
  if (1==complete) {
    // Create a new empty event.
    LADEvent* emptyev=getLADEvent(status);
    CHECK_STATUS_RET(*status, 0);
    copyLADEvent(event, emptyev);
    freeLADEvent(&emptyev);

    // Construct a combined event for output.
    event->panel   = first->signal->panel;
    event->module  = first->signal->module;
    event->element = first->signal->element;
    event->anode   = first->signal->anode;
    event->time    = first->signal->time;
    event->signal  = first->signal->signal;
    // Set PH_IDs and SRC_IDs.
    long jj;
    for (jj=0; jj<NLADSIGNALPHOTONS; jj++) {
      event->ph_id[jj] =first->signal->ph_id[jj];
      event->src_id[jj]=first->signal->src_id[jj];		
    }

    // Delete the first element from the list.
    struct LADSignalListItem* next=first->next;
    freeLADSignal(&(first->signal));
    free(first);
    first=next;

    // Search the list in order to find adjacent signals.
    float maxsignal=event->signal;

    long nanodes = 
      lad->panel[event->panel]->module[event->module]->
      element[event->element]->nanodes;
    long* anodes=malloc(nanodes/2*sizeof(long));
    CHECK_NULL_RET(anodes, *status, "memory allocation for list failed", 0);
    long ii;
    for (ii=1; ii<nanodes/2; ii++) {
      anodes[ii]=-1;
    }
    anodes[0]=event->anode;

    int new=1;
    while(1==new) {
      new=0;

      struct LADSignalListItem** item=&first;
      while (NULL!=(*item)) {

	// Check if the regarded signal in the list seems to belong to 
	// the same photon event.
	// Check for the panel, module, and element.
	if (((*item)->signal->panel==event->panel) &&
	    ((*item)->signal->module==event->module) &&
	    ((*item)->signal->element==event->element)) {
	  // Check for the anode, considering the different (bottom and top) 
	  // anode lines.
	  if ((((*item)->signal->anode>=nanodes/2)&&(event->anode>=nanodes/2)) ||
	      (((*item)->signal->anode< nanodes/2)&&(event->anode< nanodes/2))) {
	    for (ii=0; ii<nanodes/2; ii++) {
	      if (anodes[ii]==-1) break;
	      if (((*item)->signal->anode==anodes[ii]-1) || 
		  ((*item)->signal->anode==anodes[ii]+1)) {
		// Add the signal to the event.
		event->signal += (*item)->signal->signal;
		if ((*item)->signal->signal>maxsignal) {
		  maxsignal=(*item)->signal->signal;
		  event->anode = (*item)->signal->anode;
		}
		long jj;
		for (jj=0; jj<NLADSIGNALPHOTONS; jj++) {
		  long kk;
		  for (kk=0; kk<NLADEVENTPHOTONS; kk++) {
		    if ((*item)->signal->ph_id[jj]==event->ph_id[kk]) break;
		    if (0==event->ph_id[kk]) {
		      event->ph_id[kk] =(*item)->signal->ph_id[jj];
		      event->src_id[kk]=(*item)->signal->src_id[jj];		
		      break;
		    }
		  }
		}
		for (ii=0; ii<nanodes/2; ii++) {
		  if (anodes[ii]==-1) {
		    anodes[ii]=(*item)->signal->anode;
		    break;
		  }
		}
	      
		// Delete the signal entry from the list.
		next=(*item)->next;
		freeLADSignal(&((*item)->signal));
		free((*item));
		(*item)=next;

		new=1;
		break;
	      }
	    }
	    // END of loop over all anodes belonging to this signal.
	  }
	} 
	// End of check for adjacent signal.
	if (1==new) break;

	// Move to the next signal in the list.
	item=&((*item)->next);
      }
      // End of loop over all signals in the cache.
    }
    // End of searching adjacent signals.

    return(1);
  }


  if (NULL!=signal) {
    // Move to the end of the list.
    struct LADSignalListItem** item=&first;
    while (NULL!=(*item)) {
      item=&((*item)->next);
    }

    // Append the new signal to the end of the list.
    (*item)=(struct LADSignalListItem*)malloc(sizeof(struct LADSignalListItem));
    CHECK_NULL_RET((*item), *status, 
		   "memory allocation for signal list entry failed", 0);
    (*item)->signal=(LADSignal*)malloc(sizeof(LADSignal));
    CHECK_NULL_RET((*item)->signal, *status, 
		   "memory allocation for signal entry failed", 0);
    (*item)->next=NULL;
    copyLADSignal((*item)->signal, signal);
  }
  return(0);
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

  // Signal list file.
  LADSignalListFile* slf=NULL;

  // Recombined event list file.
  LADEventListFile* elf=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("ladsim");
  set_toolversion("0.14");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=ladsim_getpar(&par);
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
    
    // Determine the signal list output file.
    char signallist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.SignalList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(signallist_filename, "");
    } else {
      strcpy(signallist_filename, par.Prefix);
      strcat(signallist_filename, par.SignalList);
    }

    // Determine the event list output file.
    char eventlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.EventList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(eventlist_filename, "events.fits");
    } else {
      strcpy(eventlist_filename, par.Prefix);
      strcat(eventlist_filename, par.EventList);
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

    // Determine the appropriate detector XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_LADXMLFile(xml_filename, par.XMLFile);
    CHECK_STATUS_BREAK(status);

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
      ac->nentries=1;
      ac->entry[0] = defaultAttitudeEntry();
      ac->entry[0].time = par.TIMEZERO;
      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);
      
      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      if ((ac->entry[0].time > par.TIMEZERO) || 
	  (ac->entry[ac->nentries-1].time < par.TIMEZERO+par.Exposure)) {
	status=EXIT_FAILURE;
	char msg[MAXMSG];
	sprintf(msg, "attitude data does not cover the "
		"specified period from %lf to %lf!", 
		par.TIMEZERO, par.TIMEZERO+par.Exposure);
	SIXT_ERROR(msg);
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


    // --- Open and set up files ---

    // Open the output photon list file.
    if (strlen(photonlist_filename)>0) {
      plf=openNewPhotonListFile(photonlist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list file.
    if (strlen(impactlist_filename)>0) {
      ilf=openNewLADImpactListFile(impactlist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output raw event list file.
    if (strlen(signallist_filename)>0) {
      slf=openNewLADSignalListFile(signallist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list file for recombined events.
    elf=openNewLADEventListFile(eventlist_filename, par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    // If this is a pointing attitude, store the direction in the output
    // photon list.
    if (1==ac->nentries) {
      // Determine the telescope pointing direction and roll angle.
      Vector pointing=getTelescopeNz(ac, par.TIMEZERO, &status);
      CHECK_STATUS_BREAK(status);
    
      // Direction.
      double ra, dec;
      calculate_ra_dec(pointing, &ra, &dec);
    
      // Roll angle.
      float rollangle=getRollAngle(ac, par.TIMEZERO, &status);
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

      // Signal list file.
      if (NULL!=slf) {
	fits_update_key(slf->fptr, TDOUBLE, "RA_PNT", &ra,
			"RA of pointing direction [deg]", &status);
	fits_update_key(slf->fptr, TDOUBLE, "DEC_PNT", &dec,
			"Dec of pointing direction [deg]", &status);
	fits_update_key(slf->fptr, TFLOAT, "PA_PNT", &rollangle,
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
      if (NULL!=slf) {
	fits_update_key(slf->fptr, TSTRING, "ATTITUDE", par.Attitude,
			"attitude file", &status);
      }
      fits_update_key(elf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		      "attitude file", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Timing keywords.
    double dbuffer=0.;
    // Photon list file.
    if (NULL!=plf) {
      fits_update_key(plf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		      "reference MJD", &status);
      fits_update_key(plf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		      "time offset", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Impact list file.
    if (NULL!=ilf) {
      fits_update_key(ilf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		      "reference MJD", &status);
      fits_update_key(ilf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		      "time offset", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Signal list file.
    if (NULL!=slf) {
      fits_update_key(slf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		      "reference MJD", &status);
      fits_update_key(slf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		      "time offset", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Event list file.
    fits_update_key(elf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    fits_update_key(elf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		    "time offset", &status);
    fits_update_key(elf->fptr, TDOUBLE, "EXPOSURE", &par.Exposure,
		    "exposure time [s]", &status);
    CHECK_STATUS_BREAK(status);

    // --- End of opening files ---


    // --- Simulation Process ---

    headas_chat(3, "start simulation ...\n");

    // Loop over photon generation and processing
    // till the time of the photon exceeds the requested
    // exposure time.
    // Simulation progress status (running from 0 to 1000).
    int progress=0;
    do {

      // Photon generation.
      // Get a new photon from the generation routine.
      Photon ph;
      int isph=phgen(ac, srccat, par.TIMEZERO, par.Exposure, par.MJDREF, 
		     par.dt, lad->fov_diameter, &ph, &status);
      CHECK_STATUS_BREAK(status);
      
      // If no photon has been generated, break the loop.
      if (0==isph) break;
      
      // Check if the photon still is within the requested exposure time.
      if (ph.time>par.TIMEZERO+par.Exposure) break;

      // If requested write the photon to the output file.
      if (NULL!=plf) {
	status=addPhoton2File(plf, &ph);
	CHECK_STATUS_BREAK(status);
      }

      // Photon imaging.
      LADImpact imp;
      int isimg=ladphimg(lad, ac, &ph, &imp, &status);
      CHECK_STATUS_BREAK(status);

      // If the photon is not imaged but lost in the collimator,
      // continue with the next one.
      if (0==isimg) continue;

      // If requested write the impact to the output file.
      if (NULL!=ilf) {
	addLADImpact2File(ilf, &imp, &status);	    
	CHECK_STATUS_BREAK(status);
      }

      // Photon Detection.
      LADSignal signals[5];
      int ndet=ladphdet(lad, &imp, signals);
      assert(ndet<=5);

      // If the photon is not detected but lost,
      // continue with the next one.
      if (0==ndet) continue;

      // Write the signals to the output file.
      int ii;
      for (ii=0; ii<ndet; ii++) {
	if (NULL!=slf) {
	  addLADSignal2File(slf, &(signals[ii]), &status);
	  CHECK_STATUS_BREAK(status);
	}

	// Recombine neighboring signals to events.
	LADEvent ev;
	while (ladevrecomb(lad, &(signals[ii]), &ev, &status)>0) {
	  CHECK_STATUS_BREAK(status);

	  // Add the event to the output file.
	  addLADEvent2File(elf, &ev, &status);
	  CHECK_STATUS_BREAK(status);
	}
	CHECK_STATUS_BREAK(status);
      }
      CHECK_STATUS_BREAK(status);
      // END of loop over all signals.

      // Program progress output.
      while ((int)((ph.time-par.TIMEZERO)*1000./par.Exposure)>progress) {
	progress++;
	headas_chat(2, "\r%.1lf %%", progress*1./10.);
	fflush(NULL);
      }

    } while(1); // END of photon processing loop.
    CHECK_STATUS_BREAK(status);


    // Call recombination routine one last time to finish up.
    LADEvent ev;
    while(ladevrecomb(lad, NULL, &ev, &status)>0) {
      CHECK_STATUS_BREAK(status);
      
      // Add the event to the output file.
      addLADEvent2File(elf, &ev, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Progress output.
    headas_chat(2, "\r%.1lf %%\n", 100.);
    fflush(NULL);

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeLADEventListFile(&elf, &status);
  freeLADSignalListFile(&slf, &status);
  freeLADImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);
  freeSourceCatalog(&srccat, &status);
  freeAttitudeCatalog(&ac);
  freeLAD(&lad);

  // Release HEADAS random number generator:
  HDmtFree();

  if (status==EXIT_SUCCESS) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int ladsim_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("Prefix", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the prefix for the output files!\n", 
		   status);
    return(status);
  }
  strcpy(par->Prefix, sbuffer);
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
    SIXT_ERROR("failed reading MJDREF");
    return(status);
  } 

  status=ape_trad_query_double("TIMEZERO", &par->TIMEZERO);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading TIMEZERO");
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


