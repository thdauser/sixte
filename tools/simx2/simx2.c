#include "simx2.h"


int simx2_main() 
{
  // Program parameters.
  struct Parameters par;
  
  // Detector setup.
  GenDet* det=NULL;

  // Attitude.
  AttitudeCatalog* ac=NULL;

  // Catalog of input X-ray sources.
  XRaySourceCatalog* srccat=NULL;

  // Photon list file.
  PhotonListFile* plf=NULL;

  // Impact list file.
  ImpactListFile* ilf=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("simx2");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=simx2_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Initialize HEADAS random number generator.
    HDmtInit(par.random_seed);

    // Load the detector configuration.
    det=newGenDet(par.xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Assign the event list file to the detector data structure.
    GenDetNewEventFile(det, par.eventlist_filename, &status);
    if (EXIT_SUCCESS!=status) break;

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
      srccat = newXRaySourceCatalog(&status);
      CHECK_STATUS_BREAK(status);

      status=EXIT_FAILURE;
      headas_printf("Error: source specification on the command line not "
		    "is not implemented yet!\n");
      break;

    }
    // END of setting up the source catalog.

    // Set up photon list file.
    // TODO
    // Template for the photon list file.
    //    char photonlist_template[] = "/home/schmid/share/sixt/templates/photonlist.tpl";

    // Remove the old photon list file.
    // remove(photon_filename);

    // TODO 
    // Set up impact list file.


    // --- End of Initialization ---


    // --- Simulation Process ---

    // Photon Generation.
    phgen(det, ac, srccat, plf, par.t0, par.t0+par.exposure, &status);
    CHECK_STATUS_BREAK(status);

    // TODO Reset internal line counter of photon list file.

    // Photon Imaging.
    phimg(det, ac, plf, ilf, par.t0, par.t0+par.exposure, &status);
    CHECK_STATUS_BREAK(status);

    // TODO Reset internal line counter of impact list file.

    // Photon Detection.
    // TODO

    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);
  freeXRaySourceCatalog(&srccat);
  freeAttitudeCatalog(&ac);
  destroyGenDet(&det, &status);

  // Release HEADAS random number generator:
  HDmtFree();

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int simx2_getpar(struct Parameters* const par)
{
  int status = EXIT_SUCCESS; // Error status.

  // Get the filename of the detector XML definition file.
  if ((status = PILGetFname("xml_filename", par->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector " 
		   "XML definition file!\n", status);
    return(status);
  }

  // Get the filename of the attitude file (FITS file).
  if ((status = PILGetString("attitude_filename", par->attitude_filename))) {
    HD_ERROR_THROW("Error reading the name of the attitude file!\n", status);
    return(status);

  } else if (0==strlen(par->attitude_filename)) {
    // If the attitude filename is empty, read pointing parameters 
    // from the command line / PIL.
    if ((status = PILGetReal4("pointing_ra", &par->pointing_ra))) {
      HD_ERROR_THROW("Error reading the right ascension of the telescope "
		     "pointing direction!\n", status);
      return(status);
    }
    if ((status = PILGetReal4("pointing_dec", &par->pointing_dec))) {
      HD_ERROR_THROW("Error reading the declination of the telescope "
		     "pointing direction!\n", status);
      return(status);
    }
  }

  // Determine the name of the file that contains the input source catalog.
  // The file must have the SIMPUT format.
  if ((status = PILGetString("simput_filename", par->simput_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input source "
		   "catalog (SIMPUT)!\n", status);
    return(status);
  } 
    
  // If the source catalog filename is empty, read source parameters
  // from the command line / PIL.
  if (0==strlen(par->simput_filename)) {  
    // Determine the source type.
    if ((status = PILGetInt("src_type", &par->src_type))) {
      HD_ERROR_THROW("Error reading the source type!\n", status);
      return(status);
    } 

    // Point source.
    if (0==par->src_type) {
      if ((status = PILGetReal4("src_ra", &par->src_ra))) {
	HD_ERROR_THROW("Error reading the right ascension "
		       "of the source!\n", status);
	return(status);
      }
      if ((status = PILGetReal4("src_dec", &par->src_dec))) {
	HD_ERROR_THROW("Error reading the declination "
		       "of the source!\n", status);
	return(status);
      }
    } 
    
    // Image source.
    if (2==par->src_type) {
      if ((status = PILGetFname("src_image_filename", 
				 par->src_image_filename))) {
	HD_ERROR_THROW("Error reading the filename of the source image!\n",
		       status);
	return(status);
      }   
    }
    
    // Source flux.
    if ((status = PILGetReal4("src_flux", &par->src_flux))) {
      HD_ERROR_THROW("Error reading the source flux!\n", status);
      return(status);
    }

    // Spectrum.
    if ((status = PILGetString("src_spectrum_filename", 
			       par->src_spectrum_filename))) {
      HD_ERROR_THROW("Error reading the filename of the source spectrum!\n",
		     status);
      return(status);
    }
    // Mono-energetic source.
    if (0==strlen(par->src_spectrum_filename)) {
      if ((status = PILGetReal4("src_mono_energy", &par->src_mono_energy))) {
	HD_ERROR_THROW("Error reading the source energy "
		       "(mono-energetic source)!\n", status);
	return(status);
      }
    }
    
  }
  // END of reading the source data from the command line.

  // Get the start time of the simulation.
  if ((status = PILGetReal("t0", &par->t0))) {
    HD_ERROR_THROW("Error reading the start time t_0!\n", status);
    return(status);
  }

  // Get the exposure time for the simulation.
  if ((status = PILGetReal("exposure", &par->exposure))) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  }

  // Get the filename of the output event list file.
  if ((status = PILGetFname("eventlist_filename", par->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output " 
		   "event list file!\n", status);
    return(status);
  }

  // Get the seed for the random number generator.
  if ((status = PILGetInt("random_seed", &par->random_seed))) {
    HD_ERROR_THROW("Error reading the seed for the random "
		   "number generator!\n", status);
    return(status);
  }

  // Get the name of the FITS template directory from the 
  // environment variable.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(par->fits_templates, buffer);
  } else {
    // Could not read the environment variable.
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error reading the environment variable containing "
		   "the location of the FITS templates!\n", status);
    return(status);
  }

  return(status);
}


