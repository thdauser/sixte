#include "comadetector.h"


CoMaDetector* getCoMaDetector(struct CoMaDetectorParameters* parameters,
			       int* status)
{
  // Allocate memory for the CoMaDetector object.
  CoMaDetector* det=(CoMaDetector*)malloc(sizeof(CoMaDetector));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for CoMaDetector "
		   "object!\n", *status);
    return(det);
  }
  // Set initial values.
  det->pixels    = NULL;
  det->eventfile = NULL;

  // Call the initialization routines of the underlying data structures.
  det->pixels = getSquarePixels(&parameters->pixels, status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Create and open new event list file.
  det->eventfile = openNewCoMaEventFile(parameters->eventfile_filename,
					parameters->eventfile_template,
					status);
  if (EXIT_SUCCESS!=*status) return(det);

  return(det);
}



void freeCoMaDetector(CoMaDetector* det) 
{
  if (NULL!=det) {
    // Call the clean-up routines of the underlying data structures.
    freeSquarePixels(det->pixels);
    closeCoMaEventFile(det->eventfile);
    // Free the allocated memory.
    free(det);
  }
}



int addImpact2CoMaDetector(CoMaDetector* det, Impact* impact)
{
  int status=EXIT_SUCCESS;

  // Determine the affected pixel.
  int x, y; // Detector RAWX and RAWY.
  if (getSquarePixel(det->pixels, impact->position, &x, &y)>0) {
    CoMaEvent event = { .time   = impact->time,
			.energy = impact->energy,
			.xi     = x,
			.yi     = y };

    // Add event to event file.
    status = addCoMaEvent2File(det->eventfile, &event);
    if (EXIT_SUCCESS!=status) return(status);
  }

  return(status);
}


