#include "xmsdetector.h"


int initXMSDetector(XMSDetector* xd, struct XMSDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&xd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  status = initSquarePixels(&xd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the XMS configuration.
  // --- Currently nothing to do. ---

  // Create a new FITS event file and open it.
  status = openNewXMSEventFile(&xd->eventlist, parameters->eventlist_filename,
			       parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupXMSDetector(XMSDetector* xd)
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  cleanupSquarePixels(&xd->pixels);
  status = closeXMSEventFile(&xd->eventlist);

  return(status);
}




int addImpact2XMSDetector(XMSDetector* xd, Impact* impact)
{
  int status=EXIT_SUCCESS;


  return(status);
}


