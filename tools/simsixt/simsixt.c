#include "simsixt.h"


int simsixt_main() 
{
  // Program parameters.
  struct Parameters par;
  
  // Detector setup.
  GenDet* det=NULL;

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
  set_toolname("simsixt");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=simsixt_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Initialize HEADAS random number generator.
    HDmtInit(par.random_seed);

    // Load the detector configuration.
    det=newGenDet(par.xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set up the Attitude.
    if (strlen(par.attitude_filename)>0) {
    // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.attitude_filename,
			     par.t0, par.exposure, &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // Set up the a simple pointing attitude.

      // First allocate memory.
      ac=getAttitudeCatalog(&status);
      CHECK_STATUS_BREAK(status);

      ac->entry=(AttitudeEntry*)malloc(2*sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("memory allocation for AttitudeCatalog failed!\n", status);
	break;
      }

      // Set the values of the entries.
      ac->nentries=2;
      ac->entry[0] = defaultAttitudeEntry();
      ac->entry[1] = defaultAttitudeEntry();
      
      ac->entry[0].time = par.t0;
      ac->entry[1].time = par.t0+par.exposure;

      ac->entry[0].nz = unit_vector(par.pointing_ra*M_PI/180., 
				    par.pointing_dec*M_PI/180.);
      ac->entry[1].nz = ac->entry[0].nz;

      Vector ny = {0., 1., 0.}; // TODO
      ac->entry[0].nx = vector_product(ac->entry[0].nz, ny);
      ac->entry[1].nx = ac->entry[0].nx;
    }
    // END of setting up the attitude.

    // Set up the source catalog.
    if (strlen(par.simput_filename)>0) {
      // Load the SIMPUT X-ray source catalog.
      srccat = loadSourceCatalog(par.simput_filename, det, &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // An individual source is given on the command line.
      
      // Get an empty source catalog.
      srccat = newSourceCatalog(&status);
      CHECK_STATUS_BREAK(status);

      status=EXIT_FAILURE;
      headas_printf("Error: source specification on the command line not "
		    "is not implemented yet!\n");
      break;

    }
    // END of setting up the source catalog.


    // --- End of Initialization ---


    // --- Simulation Process ---

    headas_chat(3, "start photon generation ...\n");

    // Set up photon list file.
    // Template for the photon list file.
    char photonlist_template[MAXFILENAME];
    char photonlist_filename[MAXFILENAME];
    strcpy(photonlist_template, par.fits_templates);
    strcat(photonlist_template, "/photonlist.tpl");

    if (strlen(par.photonlist_filename)>0) {
      strcpy(photonlist_filename, par.photonlist_filename);
    } else {
      // No photon list file has been given. Therefore create a temporary 
      // file, which resides in the memory and is removed after the program 
      // terminates.
      // Open with prefix 'mem://'. 
      // TODO
      strcpy(photonlist_filename, "mem://photons.fits");
    }

    // Open the output photon list file.
    plf=openNewPhotonListFile(photonlist_filename, photonlist_template, &status);
    CHECK_STATUS_BREAK(status);
    // Set the attitude filename in the photon list (obsolete).
    char buffer[MAXMSG];
    strcpy(buffer, par.attitude_filename);
    fits_update_key(plf->fptr, TSTRING, "ATTITUDE", buffer,
		    "attitude file", &status);
    CHECK_STATUS_BREAK(status);


    // Photon Generation.
    phgen(det, ac, srccat, plf, par.t0, par.t0+par.exposure, &status);
    CHECK_STATUS_BREAK(status);


    headas_chat(3, "start photon imaging ...\n");

    // Reset internal line counter of photon list file.
    plf->row=0;

    // Set up impact list file.
    // Template for the impact list file.
    char impactlist_template[MAXFILENAME];
    char impactlist_filename[MAXFILENAME];
    strcpy(impactlist_template, par.fits_templates);
    strcat(impactlist_template, "/impactlist.tpl");

    if (strlen(par.impactlist_filename)>0) {
      strcpy(impactlist_filename, par.impactlist_filename);

      // If the old file should be removed, when it already exists,
      // the filename has to start with an exclamation mark ('!').

    } else {
      // No impact list file has been given. Therefore create a temporary 
      // file, which resides in the memory and is removed after the program 
      // terminates.
      // Open with prefix 'mem://'.
      strcpy(impactlist_filename, "mem://impacts.fits");
    }

    // Open the output impact list file.
    ilf=openNewImpactListFile(impactlist_filename, impactlist_template, &status);
    CHECK_STATUS_BREAK(status);
    // Set the attitude filename in the impact list (obsolete).
    strcpy(buffer, par.attitude_filename);
    fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", buffer,
		    "attitude file", &status);
    CHECK_STATUS_BREAK(status);


    // Photon Imaging.
    phimg(det, ac, plf, ilf, par.t0, par.t0+par.exposure, &status);
    CHECK_STATUS_BREAK(status);

    // Close the photon list file in order to save memory.
    freePhotonListFile(&plf, &status);
 

    headas_chat(3, "start photon detection ...\n");

    // Reset internal line counter of impact list file.
    ilf->row=0;

    // Set up the event list file.
    // Template for the event list file.
    char eventlist_template[MAXFILENAME];
    char eventlist_filename[MAXFILENAME];
    strcpy(eventlist_filename, par.eventlist_filename);
    strcpy(eventlist_template, par.fits_templates);
    strcat(eventlist_template, "/eventlist.tpl");

    // Open the output event list file.
    // If the old file should be removed, when it already exists,
    // the filename has to start with an exclamation mark ('!').
    elf=openNewEventListFile(eventlist_filename, eventlist_template, &status);
    CHECK_STATUS_BREAK(status);
    // Set the attitude filename in the event list (obsolete).
    strcpy(buffer, par.attitude_filename);
    fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", buffer,
		    "attitude file", &status);
    CHECK_STATUS_BREAK(status);


    // Photon Detection.
    phdetGenDet(det, ilf, elf, par.t0, par.exposure, &status);
    CHECK_STATUS_BREAK(status);


    // Close the impact list file in order to save memory.
    freeImpactListFile(&ilf, &status);

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeEventListFile(&elf, &status);
  freeImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);
  freeSourceCatalog(&srccat);
  freeAttitudeCatalog(&ac);
  destroyGenDet(&det, &status);

  // Release HEADAS random number generator:
  HDmtFree();

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int simsixt_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 


  // Get the filename of the detector XML definition file.
  status=ape_trad_query_file_name("xml_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the detector " 
		   "XML definition file!\n", status);
    return(status);
  }
  strcpy(par->xml_filename, sbuffer);
  free(sbuffer);

  // Get the filename of the attitude file (FITS file).
  status=ape_trad_query_string("attitude_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the attitude file!\n", status);
    return(status);
  } 
  strcpy(par->attitude_filename, sbuffer);
  free(sbuffer);


  if ((0==strlen(par->attitude_filename)) || 
      (0==strcmp(par->attitude_filename, " "))) {
    // If the attitude filename is empty, read pointing parameters 
    // from the command line / PIL.
    status=ape_trad_query_float("pointing_ra", &par->pointing_ra);
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error reading the right ascension of the telescope "
		     "pointing direction!\n", status);
      return(status);
    }
    status=ape_trad_query_float("pointing_dec", &par->pointing_dec);
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error reading the declination of the telescope "
		     "pointing direction!\n", status);
      return(status);
    }
  }

  // Determine the name of the file that contains the input source catalog.
  // The file must have the SIMPUT format.
  status=ape_trad_query_string("simput_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the filename of the input source "
		   "catalog (SIMPUT)!\n", status);
    return(status);
  } 
  strcpy(par->simput_filename, sbuffer);
  free(sbuffer);
    
  if ((0==strlen(par->simput_filename)) ||
      (0==strcmp(par->simput_filename, " "))) {
    // If the source catalog filename is empty, read source parameters
    // from the command line / PIL.

    // Determine the source type.
    status=ape_trad_query_int("src_type", &par->src_type);
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error reading the source type!\n", status);
      return(status);
    } 

    // Point source.
    if (0==par->src_type) {
      status=ape_trad_query_float("src_ra", &par->src_ra);
      if (EXIT_SUCCESS!=status) {
	HD_ERROR_THROW("Error reading the right ascension "
		       "of the source!\n", status);
	return(status);
      }
      status=ape_trad_query_float("src_dec", &par->src_dec);
      if (EXIT_SUCCESS!=status) {
	HD_ERROR_THROW("Error reading the declination "
		       "of the source!\n", status);
	return(status);
      }
    } 
    
    // Image source.
    if (2==par->src_type) {
      status=ape_trad_query_string("src_image_filename", &sbuffer);
      if (EXIT_SUCCESS!=status) {
	HD_ERROR_THROW("Error reading the filename of the source image!\n",
		       status);
	return(status);
      }
      strcpy(par->src_image_filename, sbuffer);
      free(sbuffer);
    }
    
    // Source flux.
    status=ape_trad_query_float("src_flux", &par->src_flux);
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error reading the source flux!\n", status);
      return(status);
    }

    // Spectrum.
    status=ape_trad_query_string("src_spectrum_filename", &sbuffer);
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error reading the filename of the source spectrum!\n",
		     status);
      return(status);
    }
    strcpy(par->src_spectrum_filename, sbuffer);
    free(sbuffer);

    // Mono-energetic source.
    if ((0==strlen(par->src_spectrum_filename)) ||
	(0==strcmp(par->src_spectrum_filename, " "))) {
      status=ape_trad_query_float("src_mono_energy", &par->src_mono_energy);
      if (EXIT_SUCCESS!=status) {
	HD_ERROR_THROW("Error reading the source energy "
		       "(mono-energetic source)!\n", status);
	return(status);
      }
    }
    
  }
  // END of reading the source data from the command line.

  // Get the start time of the simulation.
  status=ape_trad_query_double("t0", &par->t0);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the start time t_0!\n", status);
    return(status);
  }

  // Get the exposure time for the simulation.
  status=ape_trad_query_double("exposure", &par->exposure);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  }

  // Get the filename of the output photon list file.
  status=ape_trad_query_string("photonlist_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the output photon "
		   "list file!\n", status);
    return(status);
  }
  strcpy(par->photonlist_filename, sbuffer);
  free(sbuffer);

  // Get the filename of the output impact list file.
  status=ape_trad_query_string("impactlist_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the output impact "
		   "list file!\n", status);
    return(status);
  }
  strcpy(par->impactlist_filename, sbuffer);
  free(sbuffer);

  // Get the filename of the output event list file.
  status=ape_trad_query_file_name("eventlist_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the output " 
		   "event list file!\n", status);
    return(status);
  }
  strcpy(par->eventlist_filename, sbuffer);
  free(sbuffer);

  // Get the seed for the random number generator.
  status=ape_trad_query_int("random_seed", &par->random_seed);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the seed for the random "
		   "number generator!\n", status);
    return(status);
  }

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  if (NULL!=(sbuffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(par->fits_templates, sbuffer);
    // Note: the char* pointer returned by getenv should not
    // be modified nor free'd.
  } else {
    status=ape_trad_query_file_name("fits_templates", &sbuffer);
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error reading the path of the FITS "
		     "templates!\n", status);
      return(status);
    }
    strcpy(par->fits_templates, sbuffer);
    free(sbuffer);    
  }

  return(status);
}


