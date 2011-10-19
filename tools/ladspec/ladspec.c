#include "ladspec.h"


////////////////////////////////////
/** Main procedure. */
int ladspec_main() {
  // Program parameters.
  struct Parameters par; 

  // Detector setup.
  LAD* lad=NULL;

  // Input event list file.
  LADEventListFile* elf=NULL;

  // Output spectrum.
  long* spec=NULL;
  fitsfile* fptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ladspec");
  set_toolversion("0.02");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=ladspec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input pattern file.
    elf=openLADEventListFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the detector XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_LADXMLFile(xml_filename, par.XMLFile);
    CHECK_STATUS_BREAK(status);

    // Load the detector configuration.
    lad=getLADfromXML(xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Allocate memory for the output spectrum.
    headas_chat(5, "create empty spectrum with %ld channels ...\n",
		lad->rmf->NumberChannels);
    spec=(long*)malloc(lad->rmf->NumberChannels*sizeof(long));
    CHECK_NULL_BREAK(spec, status, "memory allocation for spectrum failed");

    // Initialize the spectrum with 0.
    long ii; 
    for (ii=0; ii<lad->rmf->NumberChannels; ii++) {
      spec[ii]=0;
    }

    // --- END of Initialization ---


    // --- Begin Spectrum Binning ---
    headas_chat(5, "create spectrum ...\n");

    // LOOP over all events in the FITS table.
    for (ii=0; ii<elf->nrows; ii++) {
      
      // Read the next event from the file.
      LADEvent event;
      getLADEventFromFile(elf, ii+1, &event, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the PHA channel.
      long pha=getEBOUNDSChannel(event.signal, lad->rmf);
      
      // Add the event to the spectrum.
      assert(pha-lad->rmf->FirstChannel>=0);
      assert(pha<lad->rmf->NumberChannels);      
      spec[pha-lad->rmf->FirstChannel]++;

    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

    // Store the spectrum in the output file.
    headas_chat(3, "store spectrum ...\n");

    // Create a new FITS-file (remove existing one before):
    remove(par.Spectrum);
    char buffer[MAXFILENAME];
    sprintf(buffer, "%s(%s%s)", par.Spectrum, SIXT_DATA_PATH, "/templates/ladspec.tpl");
    fits_create_file(&fptr, buffer, &status);
    CHECK_STATUS_BREAK(status);

    // Get column numbers.
    int cchannel, ccounts;
    fits_get_colnum(fptr, CASEINSEN, "CHANNEL", &cchannel, &status);
    fits_get_colnum(fptr, CASEINSEN, "COUNTS", &ccounts, &status);
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_update_key(fptr, TSTRING, "RESPFILE", "", "response file", &status);
    fits_update_key(fptr, TSTRING, "ANCRFILE", "", "ancillary response file", &status);
    fits_update_key(fptr, TLONG, "DETCHANS", &lad->rmf->NumberChannels,
		    "number of detector channels", &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all channels in the spectrum.
    for (ii=0; ii<lad->rmf->NumberChannels; ii++) {    
      long channel=ii+lad->rmf->FirstChannel;
      fits_write_col(fptr, TLONG, cchannel, ii, 1, 1, &channel, &status);
    }
    fits_write_col(fptr, TLONG, ccounts, 1, 1, lad->rmf->NumberChannels,
		   spec, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Close the files.
  freeLADEventListFile(&elf, &status);
  if (NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory.
  if (NULL!=spec) free(spec);
  freeLAD(&lad);


  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}


int ladspec_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list file");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Spectrum", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output spectrum file");
    return(status);
  } 
  strcpy(par->Spectrum, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the detector definition XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}



