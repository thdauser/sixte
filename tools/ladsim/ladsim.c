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
      
      // Determine the entrance position into a collimator hole ([m]).
      struct Point2d entrance_position;
      do {
	// Get a random position on the element.
	entrance_position.x=sixt_get_random_number()*element->xdim;
	entrance_position.y=sixt_get_random_number()*element->ydim;
      } while (!LADCollimatorOpen(entrance_position));

      // Determine the position on the detector according to the off-axis
      // angle and the orientation of the element ([m]).
      // TODO
      impact.position = entrance_position;
      
      
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
		     LADEventListFile* const elf,
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
      sqrt(2.*kB*lad->temperature*impact.position.x*1.e6/lad->efield + 
	   pow(sigma0*1.e6,2.)) /1.e6;

    // Drift time.
    double drifttime = impact.position.x / vD;

    // Determine the index of the closest anode strip.
    double y0 = impact.position.y/element->anodepitch;
    long center_anode = ((long)(y0+1)) -1;

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
    if (0>channel) {
      continue; // Photon is not detected.
    }

    // Determine the signal corresponding to the channel according 
    // to the EBOUNDS table.
    float signal = getEBOUNDSEnergy(channel, lad->rmf, 0);
    assert(signal>=0.);

    // Loop over adjacent anodes.
    long ii;
    for (ii=MAX(0,center_anode-2); 
	 ii<MIN(element->nanodes,center_anode+3); ii++) {

      LADEvent event;
      event.panel   = impact.panel;
      event.module  = impact.module;
      event.element = impact.element;
      event.anode   = ii;
      event.ph_id[0] = impact.ph_id;
      event.src_id[0] = impact.src_id;

      // Measured time.
      event.time = impact.time + drifttime;

      // Measured signal.
      double yi = (ii-center_anode)*1.0;
      event.signal = 
	signal*0.5*
	(gaussint((yi+element->anodepitch*0.5-y0)/(sigma*sqrt(2.)))-
	 gaussint((yi-element->anodepitch*0.5-y0)/(sigma*sqrt(2.))));

      // Apply thresholds.
      if (NULL!=lad->threshold_readout_lo_keV) {
	if (event.signal < *(lad->threshold_readout_lo_keV)) continue;
      }
      if (NULL!=lad->threshold_readout_up_keV) {
	if (event.signal > *(lad->threshold_readout_up_keV)) continue;
      }

      // Append the new event to the file.
      addLADEvent2File(elf, &event, status);
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

  // Event list file.
  LADEventListFile* elf=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("ladsim");
  set_toolversion("0.01");


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
    // Check the available missions, instruments, and modes.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.XMLFile);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Determine the base directory containing the XML
      // definition files.
      strcpy(xml_filename, par.data_path);
      strcat(xml_filename, "/instruments/loft/lad.xml");
    } else {
      // The XML filename has been given explicitly.
      strcpy(xml_filename, par.XMLFile);
    }
    // END of determine the XML filename.


    // Determine the photon list output file and the file template.
    char photonlist_template[MAXFILENAME];
    char photonlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.PhotonList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(photonlist_filename, par.OutputStem);
      strcat(photonlist_filename, "_photons.fits");
    } else {
      strcpy(photonlist_filename, par.PhotonList);
    }
    strcpy(photonlist_template, par.data_path);
    strcat(photonlist_template, "/templates/photonlist.tpl");

    // Determine the impact list output file and the file template.
    char impactlist_template[MAXFILENAME];
    char impactlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.ImpactList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(impactlist_filename, par.OutputStem);
      strcat(impactlist_filename, "_impacts.fits");
    } else {
      strcpy(impactlist_filename, par.ImpactList);
    }
    strcpy(impactlist_template, par.data_path);
    strcat(impactlist_template, "/templates/ladimpactlist.tpl");
    
    // Determine the event list output file and the file template.
    char eventlist_template[MAXFILENAME];
    char eventlist_filename[MAXFILENAME];
    strcpy(ucase_buffer, par.EventList);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer,"NONE")) {
      strcpy(eventlist_filename, par.OutputStem);
      strcat(eventlist_filename, "_events.fits");
    } else {
      strcpy(eventlist_filename, par.EventList);
    }
    strcpy(eventlist_template, par.data_path);
    strcat(eventlist_template, "/templates/ladeventlist.tpl");

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


    // --- End of Initialization ---


    // --- Simulation Process ---

    // Open the output photon list file.
    plf=openNewPhotonListFile(photonlist_filename, 
			      photonlist_template, 
			      &status);
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
    ilf=openNewLADImpactListFile(impactlist_filename, 
				 impactlist_template, 
				 &status);
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


    // Open the output event list file.
    elf=openNewLADEventListFile(eventlist_filename, 
				eventlist_template, 
				&status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(elf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);
    fits_update_key(elf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    dbuffer=0.;
    fits_update_key(elf->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		    "time offset", &status);
    CHECK_STATUS_BREAK(status);

    // Photon Detection.
    headas_chat(3, "start photon detection ...\n");
    ladphdet(lad, ilf, elf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);


    // Close the impact list file in order to save memory.
    freeLADImpactListFile(&ilf, &status);


    // Run the event projection.
    headas_chat(5, "start sky projection ...\n");
    // TODO
    //phproj(det, ac, elf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);


    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeLADEventListFile(&elf, &status);
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


  // Get the name of the directory containing the data
  // required for the simulations from the environment variable.
  if (NULL!=(sbuffer=getenv("SIXT_DATA_PATH"))) {
    strcpy(par->data_path, sbuffer);
    // Note: the char* pointer returned by getenv should not
    // be modified nor free'd.
  } else {
    status = EXIT_FAILURE;
    SIXT_ERROR("could not read environment variable 'SIXT_DATA_PATH'");
    return(status);
  }

  return(status);
}


