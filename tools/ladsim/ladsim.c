#include "ladsim.h"


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
  ImpactListFile* ilf=NULL;

  // Event list file.
  EventListFile* elf=NULL;

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
    strcat(impactlist_template, "/templates/impactlist.tpl");
    
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
    strcat(eventlist_template, "/templates/eventlist.tpl");

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
    // TODO
    //srccat = loadSourceCatalog(par.Simput, det, &status);
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
    // TODO
    //phgen(det, ac, srccat, plf, t0, par.Exposure, par.MJDREF, &status);
    CHECK_STATUS_BREAK(status);

    // Free the source catalog in order to save memory.
    freeSourceCatalog(&srccat, &status);
    CHECK_STATUS_BREAK(status);

    // Reset internal line counter of photon list file.
    plf->row=0;


    // Open the output impact list file.
    ilf=openNewImpactListFile(impactlist_filename, 
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
    // TODO
    //phimg(det, ac, plf, ilf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

    // Close the photon list file in order to save memory.
    freePhotonListFile(&plf, &status);
 

    // Reset internal line counter of impact list file.
    ilf->row=0;


    // Open the output event list file.
    elf=openNewEventListFile(eventlist_filename, 
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
    // TODO
    //phdetGenDet(det, ilf, elf, t0, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);


    // Close the impact list file in order to save memory.
    freeImpactListFile(&ilf, &status);


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
  freeEventListFile(&elf, &status);
  freeImpactListFile(&ilf, &status);
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


