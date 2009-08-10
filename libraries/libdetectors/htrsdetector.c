#include "htrsdetector.h"


int initHTRSDetector(HTRSDetector* xd, struct HTRSDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&xd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  //  status = initHexagonalPixels(&xd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the HTRS configuration.
  // --- Currently nothing to do. ---

  // Create a new FITS event file and open it.
  status = openNewHTRSEventFile(&xd->eventlist, parameters->eventlist_filename,
				parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupHTRSDetector(HTRSDetector* xd)
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  //  cleanupSquarePixels(&xd->pixels);
  status = closeHTRSEventFile(&xd->eventlist);

  return(status);
}

