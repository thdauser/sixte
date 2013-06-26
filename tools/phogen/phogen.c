#include "phogen.h"


int phogen_main() 
{
  // Program parameters.
  struct Parameters par;

  // Instrument setup.
  GenInst* inst=NULL;

  // Attitude.
  Attitude* ac=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat=NULL;

  // Photon list file.
  PhotonListFile* plf=NULL;
  
  // Error status.
  int status=EXIT_SUCCESS;  

  // Register HEATOOL.
  set_toolname("phogen");
  set_toolversion("0.05");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=phogen_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Determine the photon list output file.
    char photonlist_filename[MAXFILENAME];
    strcpy(photonlist_filename, par.PhotonList);

    // Determine the random number seed.
    int seed;
    if (-1!=par.Seed) {
      seed = par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed = (int)time(NULL);
    }

    // Initialize the random number generator.
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile, 
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, &status);
    CHECK_STATUS_BREAK(status);
    
    // Set up the Attitude.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
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
      ac->entry[0] = defaultAttitudeEntry();      
      ac->entry[0].time = par.TSTART;
      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);

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

    // Load the SIMPUT X-ray source catalog.
    srccat=loadSourceCatalog(par.Simput, inst->tel->arf, &status);
    CHECK_STATUS_BREAK(status);

    // --- End of Initialization ---


    // --- Photon Generation Process ---

    // Open the output photon list file.
    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    plf=openNewPhotonListFile(photonlist_filename,
			      telescop, instrume, "Normal",
			      par.MJDREF, 0.0, par.TSTART, par.TSTART+par.Exposure,
			      par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par.Attitude,
		    "attitude file", &status);

    // Start the actual photon generation (after loading required data):
    headas_chat(3, "start photon generation process ...\n");

    // Loop over photon generation, till the time of the photon exceeds 
    // the requested exposure time.
    // Simulation progress status (running from 0 to 1000).
    int progress=0;
    do {

      // Photon generation.
      Photon ph;
      int isph=phgen(ac, &srccat, 1, 
		     par.TSTART, par.TSTART+par.Exposure, 
		     par.MJDREF, par.dt, 
		     inst->tel->fov_diameter, &ph, &status);
      CHECK_STATUS_BREAK(status);

      // If no photon has been generated, break the loop.
      if (0==isph) break;

      // Check if the photon still is within the requested exposre time.
      if (ph.time>par.TSTART+par.Exposure) break;

      // Write the photon to the output file.
      status=addPhoton2File(plf, &ph);
      CHECK_STATUS_BREAK(status);

      // Program progress output.
      while ((int)((ph.time-par.TSTART)*1000./par.Exposure)>progress) {
	progress++;
	headas_chat(2, "\r%.1lf %%", progress*1./10.);
	fflush(NULL);
      }

    } while(1); // END of photon generation loop.
    CHECK_STATUS_BREAK(status);
 
    // Progress output.
    headas_chat(2, "\r%.1lf %%\n", 100.);
    fflush(NULL);

    // --- End of photon generation ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freePhotonListFile(&plf, &status);
  freeSourceCatalog(&srccat, &status);
  freeAttitude(&ac);
  destroyGenInst(&inst, &status);

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}



int phogen_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the photon list");
    return(status);
  } 
  strcpy(par->PhotonList, sbuffer);
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

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the attitude");
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the right ascension of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the declination of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_file_name("Simput", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the SIMPUT file");
    return(status);
  } 
  strcpy(par->Simput, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading MJDREF");
    return(status);
  } 

  status=ape_trad_query_double("dt", &par->dt);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading dt");
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


