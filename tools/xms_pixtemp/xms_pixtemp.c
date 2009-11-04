#include "xms_pixtemp.h"


int xms_pixtemp_main() {
  struct Parameters parameters;
  XMSEventFile eventfile;

  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("xms_pixtemp");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = xms_pixtemp_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Initialize HEADAS random number generator.
    HDmtInit(1);

    // Open the event file.
    status=openXMSEventFile(&eventfile, parameters.eventlist_filename, READWRITE);
    if (EXIT_SUCCESS!=status) break;

    // Read the EBOUNDS from the detector response file.
    struct RMF* rmf = loadRMF(parameters.rsp_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // Loop over all events in the event file.
    while((EXIT_SUCCESS==status) && (0==EventFileEOF(&eventfile.generic))) {
      
      // Read the next event from the FITS file.
      XMSEvent event;
      status=XMSEventFile_getNextRow(&eventfile, &event);
      if(EXIT_SUCCESS!=status) break;

      if ((1==event.array) && // Only events from the inner array.
	  (event.xi == parameters.pixx) && (event.xi == parameters.pixx)) {
	printf(" %lf\t%lf\n", event.time, getEnergy(event.pha, rmf));
      }

    } // End of loop over all events in the event file
    
  } while(0); // End of error handling loop


  // --- Clean Up ---

  // Close the event file.
  closeXMSEventFile(&eventfile);

  // Release HEADAS random number generator.
  HDmtFree();

  return(status);
}



int xms_pixtemp_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input file!\n", status);
  }

  else if ((status = PILGetFname("rsp_filename", parameters->rsp_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector response file!\n", status);
  }

  else if ((status = PILGetInt("pixx", &parameters->pixx))) {
    HD_ERROR_THROW("Error x-coordinate of the pixel to be analysed!\n", status);
  }

  else if ((status = PILGetInt("pixy", &parameters->pixy))) {
    HD_ERROR_THROW("Error y-coordinate of the pixel to be analysed!\n", status);
  }

  return(status);
}





