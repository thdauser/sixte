#include "htrsdetector.h"


int initHTRSDetector(HTRSDetector* hd, struct HTRSDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&hd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  status = initHexagonalPixels(&hd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the HTRS configuration.
  // --- Currently nothing to do. ---

  // Create a new FITS event file and open it.
  status = openNewHTRSEventFile(&hd->eventlist, parameters->eventlist_filename,
				parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupHTRSDetector(HTRSDetector* hd)
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  cleanupHexagonalPixels(&hd->pixels);
  status = closeHTRSEventFile(&hd->eventlist);

  return(status);
}

