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

#include "ladsim.h"


static inline LADImpact* ladphimg(const LAD* const lad,
				  Attitude* const ac,
				  Photon* const ph,
				  int* const status)
{
  assert(ph!=NULL);

  // Calculate the minimum cos-value for sources inside the FOV: 
  // (angle(x0,source) <= 1/2 * diameter)
  const double fov_min_align=cos(lad->fov_diameter/2.); 

  // Determine telescope pointing direction at the current time.
  struct Telescope telescope;
  telescope.nz=getTelescopeNz(ac, ph->time, status);
  CHECK_STATUS_RET(*status, NULL);

  // Compare the photon direction to the direction of the telescope
  // axis to check whether the photon is inside the FOV.
  Vector photon_direction=unit_vector(ph->ra, ph->dec);
  if (check_fov(&photon_direction, &telescope.nz, fov_min_align)==0) {
    // Photon is inside the FOV!
    
    // Determine telescope data like pointing direction (attitude) etc.
    // The telescope coordinate system consists of an x-, y-, and z-axis.
    getTelescopeAxes(ac, &telescope.nx, &telescope.ny, &telescope.nz, 
		     ph->time, status);
    CHECK_STATUS_RET(*status, NULL);

    // Calculate the off-axis angle ([rad]).
    double cos_theta=scalar_product(&telescope.nz, &photon_direction);

    // Apply the geometric vignetting corresponding to the projected
    // detector surface.
    double p=sixt_get_random_number(status);
    CHECK_STATUS_RET(*status, NULL);
    if (p > cos_theta) {
      // Photon is not detected, since it misses the detector.
      return(NULL);
    }

    // Apply the collimator vignetting function (if available).
    if (NULL!=lad->vignetting) {
      double theta=acos(cos_theta);
      p=sixt_get_random_number(status);
      CHECK_STATUS_RET(*status, NULL);
      if (p > get_Vignetting_Factor(lad->vignetting, ph->energy, theta, 0.)) {
	// The photon does not hit the detector at all,
	// because it is absorbed by the collimator walls.
	return(NULL);
      }
    }

    // New impact.
    LADImpact* imp=getLADImpact(status);
    CHECK_STATUS_RET(*status, NULL);    

    imp->time  =ph->time;
    imp->energy=ph->energy;
    imp->ph_id =ph->ph_id;
    imp->src_id=ph->src_id;
    
    // Determine the photon impact position on the detector:
    // Randomly select a panel, module, and element.
    imp->panel=
      (long)(sixt_get_random_number(status)*lad->npanels);
    CHECK_STATUS_RET(*status, NULL);
    imp->module=
      (long)(sixt_get_random_number(status)*lad->panel[imp->panel]->nmodules);
    CHECK_STATUS_RET(*status, NULL);
    imp->element=
      (long)(sixt_get_random_number(status)*
	     lad->panel[imp->panel]->module[imp->module]->nelements);
    CHECK_STATUS_RET(*status, NULL);
    
    // Pointer to the element.
    LADElement* element=
      lad->panel[imp->panel]->module[imp->module]->element[imp->element];
    
    // Determine the sensitive area of the element.
    float xwidth=element->xdim - 2.*element->xborder;
    float ywidth=element->ydim - 2.*element->yborder;
    
    // Determine the entrance position into a collimator hole ([m]).
    long col1, row1;
    struct Point2d entrance_position;
    do {
      // Get a random position on the sensitive area of the element.
      entrance_position.x=sixt_get_random_number(status)*xwidth;
      CHECK_STATUS_RET(*status, NULL);
      entrance_position.y=sixt_get_random_number(status)*ywidth;
      CHECK_STATUS_RET(*status, NULL);
	
      // Check if the random position lies within a hole opening.
      LADCollimatorHoleIdx(entrance_position, &col1, &row1);

    } while ((col1<0)||(row1<0));

    // If a vignetting function is given, the impact position on the 
    // detector is set equivalent to the entrance position.
    if (NULL!=lad->vignetting) {
      imp->position.x=entrance_position.x;
      imp->position.y=entrance_position.y;
    } else {
      // If no explicit vignetting function is specified,
      // check if the impact position still is inside the hole.
      // Otherwise the photon has been absorbed by the walls of the hole.
      // Make sure that it is the SAME hole as before.
      
      // Determine the position on the detector according to the off-axis
      // angle and the orientation of the element ([m]).
      
      // Determine the photon direction with respect to the telescope
      // coordinate system.
      Vector deviation;
      deviation.x=scalar_product(&telescope.nx, &photon_direction);
      deviation.y=scalar_product(&telescope.ny, &photon_direction);
      deviation.z=scalar_product(&telescope.nz, &photon_direction);
    
      // Determine the length of the vector to reach from the entrance
      // position on top of the collimator to the bottom.
      // (The collimator has a thickness of 5 mm.)
      double length=0.005/deviation.z;
    
      // Add the off-axis deviation to the entrance position.
      imp->position.x=
	entrance_position.x + deviation.x * length;
      imp->position.y=
	entrance_position.y + deviation.y * length;
      
      long col2, row2;
      LADCollimatorHoleIdx(imp->position, &col2, &row2);
      if ((col1!=col2)||(row1!=row2)) {
	freeLADImpact(&imp);
	return(NULL);
      }
    }

    // Check if the impact position still is on the active area of the SDD.
    if ((imp->position.x<0.)||(imp->position.x>xwidth)||
	(imp->position.y<0.)||(imp->position.y>ywidth)) {
      freeLADImpact(&imp);
      return(NULL);
    }

    return(imp);
  } 
  // End of FOV check.

  return(NULL);
}


/** Detection routine for the charge signals on the LAD anodes. The
    charge distribution among the neighboring anodes implemented in
    this function follows the approach of Campana et al. (2011). The
    function has to be called with impact times in chronological
    order. */ 
static inline void ladphdet(const LAD* const lad, 
			    LADImpact* const imp,
			    const int conv_with_rmf,
			    LADSignalListItem** const siglist, 
			    int* const status)
{
  // Determine the measured signal.
  float signal;

  // Check if the impact energy has to be convolved with the RMF.
  // Note that for background events drawn from an input PHA,
  // the energy must NOT be convolved with the RMF again.
  if (1==conv_with_rmf) {
    // Determine the measured detector channel (PHA channel) according 
    // to the RMF.
    // The channel is obtained from the RMF using the corresponding
    // HEAdas routine which is based on drawing a random number.
    long channel;
    returnRMFChannel(lad->rmf, imp->energy, &channel);

    // Check if the signal is really measured. If the
    // PHA channel returned by the HEAdas RMF function is '-1', 
    // the photon is not detected.
    // This can happen, if the RMF actually is an RSP, i.e., it 
    // includes ARF contributions, e.g., 
    // the detector quantum efficiency and filter transmission.
    if (channel<0) {
      headas_chat(5, "# undetected photon\n");
      return; // Photon is not detected.
    } else if (channel<lad->rmf->FirstChannel) {
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "wrong channel number: %ld", channel);
      SIXT_ERROR(msg);
      return;
    }

    // Determine the signal corresponding to the channel according 
    // to the EBOUNDS table.
    signal=getEBOUNDSEnergy(channel, lad->rmf, status);
    CHECK_STATUS_VOID(*status);

  } else {
    // The impact energy is not convolved with the RMF. Instead
    // it is directly assumed as the measured signal.
    signal=imp->energy;
  }
  assert(signal>=0.);

  // Element on the LAD.
  LADElement* element= 
    lad->panel[imp->panel]->module[imp->module]->element[imp->element];
  
  // Drift velocity.
  double vD=lad->mobility*lad->efield;
  // Drift time.
  double drifttime=imp->position.x/vD;

  // Determine the anode pitch [m].
  float anode_pitch=
    (element->ydim-2.*element->yborder)/(element->nanodes/2);

  // Determine the parameters of the charge cloud.
  // Boltzmann constant.
  const double kB=8.6173324e-5; // [eV/K]
  // Charge cloud dispersion.
  const double sigma0=20.e-6; // (according to Campana, 2011) [m]
  double sigma=
    sqrt(2.*kB*lad->temperature*imp->position.x/lad->efield+
	 pow(sigma0,2.));

  // Determine the index of the closest anode strip.
  float xwidth=element->xdim-2.*element->xborder;
  double y0;
  if (imp->position.x<0.5*xwidth) {
    y0=imp->position.y/anode_pitch;
  } else {
    y0=imp->position.y/anode_pitch+element->nanodes/2;
  }
  long anode1=((long)(y0+1.0))-1;

  // Determine the index of the second closest anode.
  long anode2=anode1;
  if ((anode1!=0) && (anode1!=element->nanodes/2) && (y0-anode1*1.0<=0.5)) {
    anode2=anode1-1;
  } else if ((anode1!=element->nanodes/2-1) && (anode1!=element->nanodes-1) &&
	     (y0-anode1*1.0>=0.5)) {
    anode2=anode1+1;
  } 

  // Determine the charge distribution according to the model 
  // of Campana et al. (2011). Distribute the charge among
  // at most 2 anodes!
  double fraction1, fraction2;
  if (anode1==anode2) {
    fraction2=0.0;
  } else if (anode1<anode2) {
    fraction2=gaussint((anode1*1.0+1.0-y0)*anode_pitch/sigma);
  } else {
    fraction2=gaussint((y0-1.0*anode1)*anode_pitch/sigma);
  }
  fraction1=1.0-fraction2;

  // Produce a signal.
  LADSignalListItem** el=siglist;

  LADSignal newsignal;
  newsignal.time     =imp->time+drifttime;    
  newsignal.panel    =imp->panel;
  newsignal.module   =imp->module;
  newsignal.element  =imp->element;
  newsignal.ph_id[0] =imp->ph_id;
  newsignal.src_id[0]=imp->src_id;
  int kk;
  for (kk=1; kk<NLADSIGNALPHOTONS; kk++) {
    newsignal.ph_id[kk] =0;
    newsignal.src_id[kk]=0;
  }

  // Main signal fraction.
  newsignal.signal=fraction1*signal;
  newsignal.anode =anode1;

  // Insert into the time-ordered list.
  while(NULL!=*el) {
    if (newsignal.time<(*el)->signal.time) break;
    el=&((*el)->next);
  }
  LADSignalListItem* newel=newLADSignalListItem(status);
  CHECK_STATUS_VOID(*status);
  copyLADSignal(&(newel->signal), &newsignal);
  newel->next=*el;
  *el=newel;

  // Secondary signal fraction.
  if (anode1!=anode2) {
    newsignal.signal=fraction2*signal;
    newsignal.anode =anode2;

    // Insert into the time-ordered list.
    while(NULL!=*el) {
      if (newsignal.time<(*el)->signal.time) break;
      el=&((*el)->next);
    }
    LADSignalListItem* newel=newLADSignalListItem(status);
    CHECK_STATUS_VOID(*status);
    copyLADSignal(&(newel->signal), &newsignal);
    newel->next=*el;
    *el=newel;
  }
  // END of loop over adjacent anodes.
}


/** Return recombined events from the detected signals. */
static inline LADEvent* ladevrecomb(const LAD* const lad,
				    LADSignal* const signal,
				    int* const status)
{
  // List of raw event signals.
  static LADSignalListItem* first=NULL;

  // Flag if an event is complete (there will be no further
  // signal contributions).
  int complete=0;
  if (NULL!=first) {
    if (NULL==signal) {
      complete=1;
    } else if (signal->time-first->signal.time>lad->coincidencetime) {
      complete=1;
    }
  }

  // Produce a new event.
  while ((1==complete) && (NULL!=first)) {
    // Create a new empty event.
    LADEvent* ev=getLADEvent(status);
    CHECK_STATUS_RET(*status, NULL);

    // Construct a combined event for output.
    ev->panel  =first->signal.panel;
    ev->module =first->signal.module;
    ev->element=first->signal.element;
    ev->anode  =first->signal.anode;
    ev->time   =first->signal.time;
    ev->signal =first->signal.signal;
    // Set PH_IDs and SRC_IDs.
    long jj;
    for (jj=0; jj<NLADSIGNALPHOTONS; jj++) {
      ev->ph_id[jj] =first->signal.ph_id[jj];
      ev->src_id[jj]=first->signal.src_id[jj];		
    }

    // Delete the first element from the list.
    LADSignalListItem* next=first->next;
    free(first);
    first=next;

    // Search the list in order to find adjacent signals.
    float maxsignal=ev->signal;

    long nanodes=
      lad->panel[ev->panel]->module[ev->module]->
      element[ev->element]->nanodes;
    long* anodes=malloc(nanodes/2*sizeof(long));
    CHECK_NULL_RET(anodes, *status, "memory allocation for list failed", NULL);
    long ii;
    for (ii=1; ii<nanodes/2; ii++) {
      anodes[ii]=-1;
    }
    anodes[0]=ev->anode;

    int new=1;
    while(1==new) {
      new=0;

      LADSignalListItem** item=&first;
      while (NULL!=(*item)) {

	// Check if the regarded signal in the list seems to belong to 
	// the same photon event.
	// Check for the panel, module, and element.
	if (((*item)->signal.panel==ev->panel) &&
	    ((*item)->signal.module==ev->module) &&
	    ((*item)->signal.element==ev->element)) {
	  // Check for the anode, considering the different (bottom and top) 
	  // anode lines.
	  if ((((*item)->signal.anode>=nanodes/2)&&(ev->anode>=nanodes/2)) ||
	      (((*item)->signal.anode< nanodes/2)&&(ev->anode< nanodes/2))) {
	    for (ii=0; ii<nanodes/2; ii++) {
	      if (anodes[ii]==-1) break;
	      if (((*item)->signal.anode==anodes[ii]-1) || 
		  ((*item)->signal.anode==anodes[ii]+1)) {
		// Add the signal to the event.
		ev->signal+=(*item)->signal.signal;
		if ((*item)->signal.signal>maxsignal) {
		  maxsignal=(*item)->signal.signal;
		  ev->anode=(*item)->signal.anode;
		}
		long jj;
		for (jj=0; jj<NLADSIGNALPHOTONS; jj++) {
		  long kk;
		  for (kk=0; kk<NLADEVENTPHOTONS; kk++) {
		    if ((*item)->signal.ph_id[jj]==ev->ph_id[kk]) break;
		    if (0==ev->ph_id[kk]) {
		      ev->ph_id[kk] =(*item)->signal.ph_id[jj];
		      ev->src_id[kk]=(*item)->signal.src_id[jj];		
		      break;
		    }
		  }
		}
		for (ii=0; ii<nanodes/2; ii++) {
		  if (anodes[ii]==-1) {
		    anodes[ii]=(*item)->signal.anode;
		    break;
		  }
		}
	      
		// Delete the signal entry from the list.
		next=(*item)->next;
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

    // Release memory.
    free(anodes);

    // Apply thresholds.
    if (NULL!=lad->threshold_readout_lo_keV) {
      if (ev->signal<=*(lad->threshold_readout_lo_keV)) {
	freeLADEvent(&ev);
	continue;
      }
    }
    if (NULL!=lad->threshold_readout_up_keV) {
      if (ev->signal>=*(lad->threshold_readout_up_keV)) {
	freeLADEvent(&ev);
	continue;
      }
    }
    return(ev);
  }


  if (NULL!=signal) {
    // Move to the end of the list.
    LADSignalListItem** item=&first;
    while (NULL!=(*item)) {
      item=&((*item)->next);
    }

    // Append the new signal to the end of the list.
    (*item)=newLADSignalListItem(status);
    CHECK_STATUS_RET(*status, NULL);
    copyLADSignal(&(*item)->signal, signal);
  }
  return(NULL);
}


int ladsim_main() 
{
  // Program parameters.
  struct Parameters par;
  
  // Detector setup.
  LAD* lad=NULL;

  // Attitude.
  Attitude* ac=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat=NULL;

  // Photon list file.
  PhotonFile* plf=NULL;

  // Impact list file.
  LADImpactFile* ilf=NULL;

  // Signal list file.
  LADSignalFile* slf=NULL;

  // Recombined event list file.
  LADEventFile* elf=NULL;

  // List of measured signals.
  LADSignalListItem* siglist=NULL;
 
  // Output file for progress status.
  FILE* progressfile=NULL;

  // Artificial flat ARF for the generation of background events.
  struct ARF* bkgarf=NULL;

  // Background source.
  SimputSrc* bkgsrc=NULL;

  // GTI collection.
  GTI* gti=NULL;

  // Number of simulated background events.
  long nbkgevts=0;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("ladsim");
  set_toolversion("0.31");


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
    strcpy(ucase_buffer, par.EvtFile);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(eventlist_filename, par.Prefix);
      strcat(eventlist_filename, "evt.fits");
    } else {
      strcpy(eventlist_filename, par.Prefix);
      strcat(eventlist_filename, par.EvtFile);
    }

    // Initialize the random number generator.
    sixt_init_rng(getSeed(par.Seed), &status);
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

    // Determine the appropriate detector XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_LADXMLFile(xml_filename, par.XMLFile);
    CHECK_STATUS_BREAK(status);

    // Load the detector configuration.
    lad=getLADfromXML(xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set up the background ARF if necessary.
    if (NULL!=lad->bkgctlg) {
      // Get an empty ARF.
      bkgarf=getARF(&status);
      CHECK_STATUS_BREAK(status);

      // Allocate memory. Use a logarithmic energy grid for this ARF.
      // So the number of required bins is determined as 
      // log(Emax/Emin)/log(factor)
      // For Emin=0.9, Emax=1000., and factor=1.01 a total number of
      // 704 bins is required.
      const long narfbins=704;
      const float arfEmin=0.9; // lower boundary for the energy grid.
      const float arffactor=1.01;
      bkgarf->LowEnergy=(float*)malloc(narfbins*sizeof(float));
      CHECK_NULL_BREAK(bkgarf->LowEnergy, status, 
		       "memory allocation for background ARF failed");
      bkgarf->HighEnergy=(float*)malloc(narfbins*sizeof(float));
      CHECK_NULL_BREAK(bkgarf->HighEnergy, status, 
		       "memory allocation for background ARF failed");
      bkgarf->EffArea=(float*)malloc(narfbins*sizeof(float));
      CHECK_NULL_BREAK(bkgarf->EffArea, status, 
		       "memory allocation for background ARF failed");
      bkgarf->NumberEnergyBins=narfbins;

      // Set up a simple ARF.
      float sensitive_area=
	lad->npanels*lad->panel[0]->nmodules*
	lad->panel[0]->module[0]->nelements*
	(lad->panel[0]->module[0]->element[0]->xdim-
	 2.*lad->panel[0]->module[0]->element[0]->xborder)*
	(lad->panel[0]->module[0]->element[0]->ydim-
	 2.*lad->panel[0]->module[0]->element[0]->yborder)*
	1.e4; // Total sensitive area [cm^2].
      long ii;
      for (ii=0; ii<narfbins; ii++) {
	bkgarf->LowEnergy[ii] =arfEmin*pow(arffactor,ii); // [keV]
	bkgarf->HighEnergy[ii]=arfEmin*pow(arffactor,ii+1);
	bkgarf->EffArea[ii]   =sensitive_area;
      }
      headas_chat(5, "background ARF with %.2f cm^2\n", sensitive_area);

      // Set reference to background ARF for SIMPUT library.
      setSimputARF(lad->bkgctlg, bkgarf);

      // Set reference to the background source.
      bkgsrc=loadSimputSrc(lad->bkgctlg, 1, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Set up the Attitude.
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a pointing attitude.
      ac=getPointingAttitude(par.MJDREF, par.TSTART, par.TSTART+par.Exposure,
			     par.RA*M_PI/180., par.Dec*M_PI/180., &status);
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

    // Load the SIMPUT X-ray source catalog.
    srccat=loadSourceCatalog(par.Simput, lad->arf, &status);
    CHECK_STATUS_BREAK(status);


#ifdef LAD_OAR
    // Determine the Open Area Ratio of the Collimator on the LAD.
    long kk;
    long try=0, pass=0;
    struct Point2d position;
    float xwidth=
      lad->panel[0]->module[0]->element[0]->xdim-
      lad->panel[0]->module[0]->element[0]->xborder*2.;
    float ywidth=
      lad->panel[0]->module[0]->element[0]->ydim-
      lad->panel[0]->module[0]->element[0]->yborder*2.;
    for (kk=0; kk<1000000; kk++) {
      long col, row;
      do {
	try++;
	// Get a random position on the element.
	position.x=sixt_get_random_number(&status)*xwidth;
	CHECK_STATUS_BREAK(status);
	position.y=sixt_get_random_number(&status)*ywidth;
	CHECK_STATUS_BREAK(status);
	
	// Determine the indices of the respective hole.
	LADCollimatorHoleIdx(position, &col, &row);

      } while ((col<0)||(row<0));
      pass++;
    }
    headas_chat(1, "### LAD Open Area Ratio: %lf ###\n", pass*1./try);
#endif

    // --- End of Initialization ---


    // --- Open and set up files ---

    // Open the output photon list file.
    if (strlen(photonlist_filename)>0) {
      plf=openNewPhotonFile(photonlist_filename, 
			    "LOFT", "LAD", "Normal",
			    lad->arf_filename, lad->rmf_filename,
			    par.MJDREF, 0.0, 
			    par.TSTART, par.TSTART+par.Exposure,
			    par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output impact list file.
    if (strlen(impactlist_filename)>0) {
      ilf=openNewLADImpactFile(impactlist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output raw event list file.
    if (strlen(signallist_filename)>0) {
      slf=openNewLADSignalFile(signallist_filename, par.clobber, &status);
      CHECK_STATUS_BREAK(status);
    }

    // Open the output event list file for recombined events.
    elf=openNewLADEventFile(eventlist_filename, 
			    lad->arf_filename, lad->rmf_filename,
			    par.MJDREF, 0.0,
			    par.TSTART, par.TSTART+par.Exposure,
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

    // Store the exposure time.
    fits_update_key(elf->fptr, TDOUBLE, "EXPOSURE", &par.Exposure,
		    "exposure time [s]", &status);
    CHECK_STATUS_BREAK(status);

    // --- End of opening files ---


    // --- Simulation Process ---

    headas_chat(3, "start simulation ...\n");

    // Simulation progress status (running from 0 to 100).
    int progress=0;
    if (NULL==progressfile) {
      headas_chat(2, "\r%.0lf %%", 0.);
      fflush(NULL);
    } else {
      rewind(progressfile);
      fprintf(progressfile, "%.2lf", 0.);
      fflush(progressfile);	
    }

    // Loop over photon generation and processing
    // till the time of the photon exceeds the requested
    // exposure time.
    double t_next_bkg=par.TSTART;
    do {
      // Get a new photon from the generation routine.
      Photon ph;
      int isph=phgen(ac, &srccat, 1, par.TSTART, par.TSTART+par.Exposure, 
		     par.MJDREF, par.dt, lad->fov_diameter, &ph, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the photon still is within the requested exposure time.
      if ((0!=isph)&&(ph.time>par.TSTART+par.Exposure)) {
	isph=0;
      }
      
      // If a photon has been generated within the regarded time interval.
      LADImpact* imp=NULL;
      if (0!=isph) {
	// If requested, write the photon to the output file.
	if (NULL!=plf) {
	  status=addPhoton2File(plf, &ph);
	  CHECK_STATUS_BREAK(status);
	}

	// Photon transmission through / absorption in the collimator.
	imp=ladphimg(lad, ac, &ph, &status);
	CHECK_STATUS_BREAK(status);

	// If the photon is not transmitted but lost in the collimator,
	// continue with the next one.
	if (NULL==imp) continue;
    
	// If requested, write the impact to the output file.
	if (NULL!=ilf) {
	  addLADImpact2File(ilf, imp, &status);	    
	  CHECK_STATUS_BREAK(status);
	}
      } else {
	// No photon has been produced within the regarded time interval.
	// So we have to finalize the simulation, i.e. get all remaining
	// events from the internal caches.
	imp=getLADImpact(&status);
	CHECK_STATUS_BREAK(status);
	imp->time=par.TSTART+par.Exposure;
	imp->energy=-1.;
      }

      while(imp!=NULL) {

	// Insert background events.
	double readouttime;
	if ((NULL!=lad->bkgctlg) && (t_next_bkg<imp->time)) {
	  LADImpact* bkgimp=getLADImpact(&status);
	  CHECK_STATUS_BREAK(status);
	  
	  // Time.
	  bkgimp->time=t_next_bkg;
      
	  // Energy.
	  double ra, dec;
	  getSimputPhotonEnergyCoord(lad->bkgctlg, bkgsrc,
				     bkgimp->time, par.MJDREF,
				     &bkgimp->energy, &ra, &dec,
				     &status);
	  CHECK_STATUS_BREAK(status);

	  // Position.
	  bkgimp->panel=
	    (long)(sixt_get_random_number(&status)*lad->npanels);
	  CHECK_STATUS_BREAK(status);
	  bkgimp->module=
	    (long)(sixt_get_random_number(&status)*
		   lad->panel[bkgimp->panel]->nmodules);
	  CHECK_STATUS_BREAK(status);
	  bkgimp->element=
	    (long)(sixt_get_random_number(&status)*
		   lad->panel[bkgimp->panel]->
		   module[bkgimp->module]->nelements);
	  CHECK_STATUS_BREAK(status);
	  // Pointer to the element.
	  LADElement* ladelement=
	    lad->panel[bkgimp->panel]->module[bkgimp->module]->
	    element[bkgimp->element];
	  // Determine the sensitive area of the element.
	  float xwidth=ladelement->xdim-2.0*ladelement->xborder;
	  float ywidth=ladelement->ydim-2.0*ladelement->yborder;
	  // Get a random position on the element.
	  bkgimp->position.x=sixt_get_random_number(&status)*xwidth;
	  CHECK_STATUS_BREAK(status);
	  bkgimp->position.y=sixt_get_random_number(&status)*ywidth;
	  CHECK_STATUS_BREAK(status);
	  
	  // Photon and source ID.
	  bkgimp->ph_id =0;
	  bkgimp->src_id=0;
	  
	  // Insert the background impact into the time-ordered cache.
	  ladphdet(lad, bkgimp, 0, &siglist, &status);
	  CHECK_STATUS_BREAK(status);
	  readouttime=bkgimp->time;

	  // Release memory.
	  freeLADImpact(&bkgimp);

	  // Count the number of background events.
	  nbkgevts++;

	  // Determine the time of the next background event.
	  int failed=
	    getSimputPhotonTime(lad->bkgctlg, bkgsrc, t_next_bkg,
				par.MJDREF, &t_next_bkg, &status);
	  CHECK_STATUS_BREAK(status);
	  if (0!=failed) {
	    SIXT_ERROR("failed producing background event");
	    status=EXIT_FAILURE;
	    break;
	  }
	} else {
	  // Insert the foreground event.
	  if (imp->energy>=0.) {
	    ladphdet(lad, imp, 1, &siglist, &status);
	    CHECK_STATUS_BREAK(status);
	  }
	  readouttime=imp->time;
	  freeLADImpact(&imp);
	}

	// Determine the signals at the individual anodes.
	while (NULL!=siglist) {
	  
	  LADSignal* signal=&siglist->signal;

	  // Go through the cache of detected signals and read out
	  // all which happend before imp->time-tDmax. We need to
	  // do this, because due to the different drift times, the time
	  // order of the impacts is shuffled.
	  if (signal->time>readouttime) {
	    break;
	  }

	  // Apply the ASIC dead time.
	  // Element on the LAD.
	  LADElement* element=
	    lad->panel[signal->panel]->module[signal->module]->
	    element[signal->element];

	  // Determine the ASIC and the pin of the ASIC the
	  // anode is attached to.
	  int asic=(int)(signal->anode/lad->asic_channels);
	  int pin =      signal->anode%lad->asic_channels;
	  
	  // Check if the pin is at the border of the ASIC and whether
	  // the neighboring ASIC has to be read out, too.
	  int asic2=-1;
	  if ((0==pin) && (asic!=0) && (asic!=element->nasics/2)) {
	    asic2=asic-1;
	  } else if ((pin>lad->asic_channels-2) &&
		     (asic!=element->nasics/2-1) &&
		     (asic!=element->nasics-1)) {
	    asic2=asic+1;
	  }

	  // Check if the event happens after the coincidence time, but
	  // during the dead time.
	  if (element->asic_readout_time[asic]>0.) {
	    if ((signal->time-element->asic_readout_time[asic]>
		 lad->coincidencetime) &&
		(signal->time-element->asic_readout_time[asic]<
		 lad->coincidencetime+element->asic_deadtime[asic])) {
	      
	      // Delete the element from the buffered list.
	      LADSignalListItem* next=siglist->next;
	      free(siglist);
	      siglist=next;
	      continue;
	    }
	  }
	  if (asic2>=0) {
	    if (element->asic_readout_time[asic2]>0.) {
	      if ((signal->time-element->asic_readout_time[asic2]>
		   lad->coincidencetime) &&
		  (signal->time-element->asic_readout_time[asic2]<
		   lad->coincidencetime+element->asic_deadtime[asic2])) {

		// Delete the element from the buffered list.
		LADSignalListItem* next=siglist->next;
		free(siglist);
		siglist=next;
		continue;
	      }
	    }
	  }
    
	  // Set the time of this ASIC readout.
	  if ((signal->time-element->asic_readout_time[asic]>
	       lad->coincidencetime) ||
	      (element->asic_readout_time[asic]==0.)) {
	    element->asic_readout_time[asic]=signal->time;
	    if (asic2>=0) {
	      element->asic_readout_time[asic2]=signal->time;
	    }
	  }
    
	  // Write the detected signal to the output file.
	  if (NULL!=slf) {
	    addLADSignal2File(slf, signal, &status);
	    CHECK_STATUS_BREAK(status);
	  }

	  // Recombine neighboring signals to events.
	  LADEvent* ev;
	  while ((ev=ladevrecomb(lad, &siglist->signal, &status))) {
	    CHECK_STATUS_BREAK(status);

	    // Add the event to the output file.
	    addLADEvent2File(elf, ev, &status);
	    CHECK_STATUS_BREAK(status);
	  
	    freeLADEvent(&ev);
	  }
	  CHECK_STATUS_BREAK(status);
	  // END of loop over all events.

	  // Move to the next entry.
	  LADSignalListItem* next=siglist->next;
	  free(siglist);
	  siglist=next;
	}
	CHECK_STATUS_BREAK(status);
	// END of loop over all signals.
      }
      CHECK_STATUS_BREAK(status);
      // END of loop while there is an impact.

      if (0==isph) break;
	
      // Program progress output.
      while ((int)((ph.time-par.TSTART)*100./par.Exposure)>progress) {
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

    } while(1); // END of photon processing loop.
    CHECK_STATUS_BREAK(status);

    // Make sure that the signal list has been processed until
    // the end of the simulated interval.
    if (siglist!=NULL) {
      assert(siglist->signal.time>par.TSTART+par.Exposure);
    }

    // Store the GTI extension in the event file.
    gti=newGTI(&status);
    CHECK_STATUS_BREAK(status);
    gti->mjdref=par.MJDREF;
    gti->timezero=0.0;
    appendGTI(gti, par.TSTART, par.TSTART+par.Exposure, &status);
    CHECK_STATUS_BREAK(status);
    saveGTIExt(elf->fptr, "STDGTI", gti, &status);
    CHECK_STATUS_BREAK(status);

    // Progress output.
    if (NULL==progressfile) {
      headas_chat(2, "\r%.0lf %%\n", 100.);
      fflush(NULL);
    } else {
      rewind(progressfile);
      fprintf(progressfile, "%.2lf", 1.);
      fflush(progressfile);	
    }

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeLADEventFile(&elf, &status);
  freeLADSignalFile(&slf, &status);
  freeLADImpactFile(&ilf, &status);
  freePhotonFile(&plf, &status);
  freeSourceCatalog(&srccat, &status);
  freeAttitude(&ac);
  freeLADSignalList(&siglist);
  freeLAD(&lad, &status);
  freeSimputSrc(&bkgsrc);
  freeARF(bkgarf);
  freeGTI(&gti);
  bkgarf=NULL;

  if (NULL!=progressfile) {
    fclose(progressfile);
    progressfile=NULL;
  }

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) {
    headas_chat(5, "number of background events: %ld\n", nbkgevts);
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
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

  status=ape_trad_query_string("RawData", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the raw event list!\n", 
		   status);
    return(status);
  } 
  strcpy(par->SignalList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the event list!\n", status);
    return(status);
  } 
  strcpy(par->EvtFile, sbuffer);
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

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the seed for the random number generator");
    return(status);
  }

  status=ape_trad_query_string("ProgressFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the progress status file!\n", 
		   status);
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


