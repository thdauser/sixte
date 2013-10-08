#include "evpat.h"


int evpat_main() 
{
  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Input event list file.
  EventFile* elf=NULL;

  // Output pattern file.
  PatternFile* plf=NULL;

  // Instrument.
  GenInst* inst=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("evpat");
  set_toolversion("0.04");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile, 
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the event list file name.
    char eventlist_filename[MAXFILENAME];
    strcpy(eventlist_filename, par.EventList);

    // Determine the pattern output file.
    char pattern_filename[MAXFILENAME];
    strcpy(pattern_filename, par.PatternList);


    headas_chat(3, "start pattern recombination ...\n");

    // Open the input event file.
    elf=openEventFile(eventlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read the timing keywords.
    char comment[MAXMSG];
    double mjdref, timezero, tstart, tstop;
    fits_read_key(elf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output pattern file.
    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    plf=openNewPatternFile(pattern_filename, 
			   telescop, instrume, "Normal",
			   inst->tel->arf_filename,
			   inst->det->rmf_filename,
			   mjdref, timezero, tstart, tstop,			   
			   inst->det->pixgrid->xwidth,
			   inst->det->pixgrid->ywidth,
			   par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Pattern recombination.
    phpat(inst->det, elf, plf, par.SkipInvalids, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&elf, &status);
  destroyPatternFile(&plf, &status);
 
  destroyGenInst(&inst, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the input event list");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output pattern list");
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
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

  status=ape_trad_query_bool("SkipInvalids", &par->SkipInvalids);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the SkipInvalids parameter");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}


