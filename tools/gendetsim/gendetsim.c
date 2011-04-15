#include "gendetsim.h"


////////////////////////////////////
/** Main procedure. */
int gendetsim_main() {

  // Containing all programm parameters read by PIL
  struct Parameters par; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;
  // Input impact list.
  ImpactListFile* ilf=NULL;

  // Output event list file.
  EventListFile* elf=NULL;

  int status=EXIT_SUCCESS; // Error status.

  // Register HEATOOL:
  set_toolname("gendetsim");
  set_toolversion("0.01");

  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

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


    // Determine the impact list output file and the file template.
    char impactlist_filename[MAXFILENAME];
    strcpy(impactlist_filename, par.ImpactList);
    
    // Determine the event list output file and the file template.
    char eventlist_template[MAXFILENAME];
    char eventlist_filename[MAXFILENAME];
    strcpy(eventlist_filename, par.EventList);
    strcpy(eventlist_template, par.fits_templates);
    strcat(eventlist_template, "/eventlist.tpl");

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

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(3, "start detection process ...\n");

    // Open the FITS file with the input impact list:
    ilf=openImpactListFile(impactlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output event file.
    elf=openNewEventListFile(eventlist_filename, eventlist_template, &status);
    CHECK_STATUS_BREAK(status);

    // Photon detection.
    phdetGenDet(det, ilf, elf, par.TIMEZERO, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Destroy the detector data structure.
  destroyGenDet(&det, &status);

  // Close the event list FITS file.
  freeEventListFile(&elf, &status);

  // Close the impact list FITS file.
  freeImpactListFile(&ilf, &status);

  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("ImpactList", &sbuffer);
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


