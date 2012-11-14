#include "gendetsim.h"


////////////////////////////////////
/** Main procedure. */
int gendetsim_main() {

  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Instrument data structure (containing the pixel array, its width, ...).
  GenInst* inst=NULL;

  // Input impact list.
  ImpactListFile* ilf=NULL;

  // Output event list file.
  EventListFile* elf=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("gendetsim");
  set_toolversion("0.02");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile,
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Use the background if available.
    setGenInstIgnoreBkg(inst, 0);

    // Set the start time for the simulation.
    setGenInstStartTime(inst, par.TIMEZERO);
    
    // Determine the impact list file.
    char impactlist_filename[MAXFILENAME];
    strcpy(impactlist_filename, par.ImpactList);
    
    // Determine the event list output file.
    char eventlist_filename[MAXFILENAME];
    strcpy(eventlist_filename, par.EventList);

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

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(3, "start detection process ...\n");

    // Open the FITS file with the input impact list:
    ilf=openImpactListFile(impactlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output event file.
    elf=openNewEventListFile(eventlist_filename, par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    if (NULL!=inst->tel->telescope) {
      fits_update_key(elf->fptr, TSTRING, "TELESCOP", inst->tel->telescope,
		      "telescope name", &status);
      CHECK_STATUS_BREAK(status);
    }
    fits_update_key(elf->fptr, TDOUBLE, "MJDREF", &par.MJDREF,
		    "reference MJD", &status);
    fits_update_key(elf->fptr, TDOUBLE, "TIMEZERO", &par.TIMEZERO,
		    "time offset", &status);
    CHECK_STATUS_BREAK(status);

    // Store the number of simulated input photons in the FITS header
    // of the output event file.
    fits_update_key(elf->fptr, TLONG, "NPHOTONS", 
		    &ilf->nrows, "number of input photons", 
		    &status);
    CHECK_STATUS_BREAK(status);

    // Define the event list file as output file.
    setGenInstEventListFile(inst, elf);

    // Loop over all impacts in the FITS file.
    while (ilf->row<ilf->nrows) {

      Impact impact;
      getNextImpactFromFile(ilf, &impact, &status);
      CHECK_STATUS_BREAK(status);

      // Check whether we are still within the requested time interval.
      if (impact.time < par.TIMEZERO) continue;
      if (impact.time > par.TIMEZERO+par.Exposure) break;

      // Photon detection.
      phdetGenDet(inst->det, &impact, par.TIMEZERO+par.Exposure, &status);
      CHECK_STATUS_BREAK(status);

    };
    CHECK_STATUS_BREAK(status);
    // End of loop over all impacts in the input file.

    // Finalize the photon detection.
    phdetGenDet(inst->det, NULL, par.TIMEZERO+par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Destroy the GenInst data structure.
  destroyGenInst(&inst, &status);

  // Close the event list FITS file.
  freeEventListFile(&elf, &status);

  // Close the impact list FITS file.
  freeImpactListFile(&ilf, &status);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully\n\n");

  return(status);
}



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

  return(status);
}


