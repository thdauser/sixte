#include "evpat.h"


int evpat_main() 
{
  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Input event list file.
  EventListFile* elf=NULL;

  // Output pattern file.
  PatternFile* plf=NULL;

  // Detector.
  GenDet* det=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("evpat");
  set_toolversion("0.02");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;

    // Determine the appropriate detector XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile, 
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Determine the event list file name.
    char eventlist_filename[MAXFILENAME];
    strcpy(eventlist_filename, par.EventList);

    // Determine the pattern output file.
    char pattern_filename[MAXFILENAME];
    strcpy(pattern_filename, par.PatternList);

    // Load the detector configuration.
    det=newGenDet(xml_filename, 1, &status);
    CHECK_STATUS_BREAK(status);


    headas_chat(3, "start pattern recombination ...\n");

    // Open the input event file.
    elf=openEventListFile(eventlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output pattern file.
    plf=openNewPatternFile(pattern_filename, par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set header keywords.
    if (NULL!=det->telescope) {
      fits_update_key(plf->fptr, TSTRING, "TELESCOP", det->telescope,
		      "telescope name", &status);
      CHECK_STATUS_BREAK(status);
    }

    // Copy FITS header keywords.
    char comment[MAXMSG];
    double mjdref=0.;
    fits_read_key(elf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    fits_update_key(plf->fptr, TDOUBLE, "MJDREF", &mjdref, 
		    "reference MJD", &status);
    CHECK_STATUS_BREAK(status);
    double timezero=0.;
    fits_read_key(elf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status);
    fits_update_key(plf->fptr, TDOUBLE, "TIMEZERO", &timezero, 
		    "time offset", &status);
    CHECK_STATUS_BREAK(status);

    // Pattern recombination.
    phpat(det, elf, plf, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventListFile(&elf, &status);
  destroyPatternFile(&plf, &status);
 
  destroyGenDet(&det, &status);

  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the input event list!\n", status);
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the output pattern list!\n", status);
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
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

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


