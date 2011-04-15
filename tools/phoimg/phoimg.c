#include "phoimg.h"


////////////////////////////////////
/** Main procedure. */
int phoimg_main() {
  struct Parameters par;

  AttitudeCatalog* ac=NULL;

  PhotonListFile* plf=NULL;
  ImpactListFile* ilf=NULL;

  // Detector data structure including telescope information like the PSF,
  // vignetting function, focal length, and FOV diameter.
  GenDet* det=NULL;

  int status=EXIT_SUCCESS; // Error status


  // Register HEATOOL:
  set_toolname("phoimg");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=phoimg_getpar(&par))) break;

    headas_chat(3, "initialize ...\n");

    // Determine the appropriate detector XML definition file.
    char xml_filename[MAXFILENAME];
    // Convert the user input to capital letters.
    strtoupper(par.Mission);
    strtoupper(par.Instrument);
    strtoupper(par.Mode);
    // Check the available missions, instruments, and modes.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.XMLFile);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Determine the base directory containing the XML
      // definition files.
      strcpy(xml_filename, par.xml_path);

      // Determine the XML filename according to the selected
      // mission, instrument, and mode.
      if (0==strcmp(par.Mission, "SRG")) {
	strcat(xml_filename, "/srg");
	if (0==strcmp(par.Instrument, "EROSITA")) {
	  strcat(xml_filename, "/erosita.xml");
	} else {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("selected instrument is not supported");
	  break;
	}

      } else if (0==strcmp(par.Mission, "IXO")) {
	strcat(xml_filename, "/ixo");
	if (0==strcmp(par.Instrument, "WFI")) {
	  strcat(xml_filename, "/wfi");
	  if (0==strcmp(par.Instrument, "FULLFRAME")) {
	    strcat(xml_filename, "/fullframe.xml");
	  } else {
	    status=EXIT_FAILURE;
	    SIXT_ERROR("selected mode is not supported");
	    break;
	  }
	} else {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("selected instrument is not supported");
	  break;
	}

      } else if (0==strcmp(par.Mission, "GRAVITAS")) {
	strcat(xml_filename, "/gravitas");
	if (0==strcmp(par.Instrument, "HIFI")) {
	  strcat(xml_filename, "/hifi.xml");
	} else {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("selected instrument is not supported");
	  break;
	}

      } else {
	status=EXIT_FAILURE;
	SIXT_ERROR("selected mission is not supported");
	break;
      }
	    
    } else {
      // The XML filename has been given explicitly.
      strcpy(xml_filename, par.XMLFile);
    }
    // END of determine the XML filename.


    // Determine the impact list input file and the file template.
    char photonlist_filename[MAXFILENAME];
    strcpy(photonlist_filename, par.PhotonList);

    // Determine the impact list output file and the file template.
    char impactlist_template[MAXFILENAME];
    char impactlist_filename[MAXFILENAME];
    strcpy(impactlist_filename, par.ImpactList);
    strcpy(impactlist_template, par.fits_templates);
    strcat(impactlist_template, "/impactlist.tpl");

    // Determine the random number seed.
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
    det=newGenDet(xml_filename, &status);
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
      
      ac->entry[0].time = par.TIMEZERO;
      ac->entry[1].time = par.TIMEZERO+par.Exposure;

      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);
      ac->entry[1].nz = ac->entry[0].nz;

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);
      ac->entry[1].nx = ac->entry[0].nx;

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, par.TIMEZERO, par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of setting up the attitude.
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(3, "start imaging process ...\n");

    // Open the input photon list file.
    plf=openPhotonListFile(photonlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output impact list file.
    ilf=openNewImpactListFile(impactlist_filename, impactlist_template, &status);
    CHECK_STATUS_BREAK(status);
    // Set the attitude filename in the impact list (obsolete).
    char buffer[MAXMSG];
    strcpy(buffer, par.Attitude);
    fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", buffer,
		    "attitude file", &status);
    CHECK_STATUS_BREAK(status);

    phimg(det, ac, plf, ilf, par.TIMEZERO, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

    // --- END of imaging process ---

  } while(0); // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the FITS files.
  freeImpactListFile(&ilf, &status);
  freePhotonListFile(&plf, &status);

  freeAttitudeCatalog(&ac);
  destroyGenDet(&det, &status);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");

  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int phoimg_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("PhotonList", &sbuffer);
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

  status=ape_trad_query_string("Mission", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the mission!\n", status);
    return(status);
  } 
  strcpy(par->Mission, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Instrument", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the instrument!\n", status);
    return(status);
  } 
  strcpy(par->Instrument, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mode", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the instrument mode!\n", status);
    return(status);
  } 
  strcpy(par->Mode, sbuffer);
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
    HD_ERROR_THROW("Error reading the right ascension of the telescope pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the declination of the telescope pointing!\n", status);
    return(status);
  } 

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


  // Get the name of the FITS template directory
  // from the environment variable.
  if (NULL!=(sbuffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(par->fits_templates, sbuffer);
    // Note: the char* pointer returned by getenv should not
    // be modified nor free'd.
  } else {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error reading the environment variable 'SIXT_FITS_TEMPLATES'!\n", 
		   status);
    return(status);
  }

  // Get the name of the directory containing the detector
  // XML definition files from the environment variable.
  if (NULL!=(sbuffer=getenv("SIXT_XML_PATH"))) {
    strcpy(par->xml_path, sbuffer);
    // Note: the char* pointer returned by getenv should not
    // be modified nor free'd.
  } else {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error reading the environment variable 'SIXT_XML_PATH'!\n", 
		   status);
    return(status);
  }

  return(status);
}



