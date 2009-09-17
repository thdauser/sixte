#include "uniformdetectorbackground.h"


int initUniformDetectorBackground(UniformDetectorBackground* background, 
				  struct UniformDetectorBackgroundParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Check if a detector background spectrum was specified.
  if (strlen(parameters->spectrum_filename) > 0) {
    // Load the background spectrum from the given PHA file.
    status = loadSpectrum(&background->spectrum, parameters->spectrum_filename);
    if (EXIT_SUCCESS!=status) return(status);
    // Set the background event rate.
    background->rate = parameters->rate;

  } else { // No background spectrum specified.
    // Set the background event rate to the default value of 0.
    background->rate = 0.;
    headas_chat(1, "Warning: no detector background spectrum specified!\n");
  }

  // Set the nextImpact (representing the next detector background event) 
  // to the default value.
  Impact emptyImpact = { .energy=0., .time=HUGE, .position={.x=0., .y=0.} };
  background->nextImpact = emptyImpact;

  return(status);
}



int cleanupUniformDetectorBackground(UniformDetectorBackground* background)
{
  // Release allocated memory.
  cleanupSpectrum(&background->spectrum);

  return(EXIT_SUCCESS);
}



int createUniformDetectorBackgroundImpact(UniformDetectorBackground* background, 
					  SquarePixels* pixels, struct RMF* rmf)
{
  // If no background is used, return without doing anything.
  if (0.==background->rate) return(EXIT_SUCCESS);

  // Create a new impact representing a detector background event.
  
  // Determine the position of the event on the CCD.
  background->nextImpact.position.x = 
    (2.*get_random_number() -1.) * pixels->xoffset * pixels->xpixelwidth;
  background->nextImpact.position.y = 
    (2.*get_random_number() -1.) * pixels->yoffset * pixels->ypixelwidth;

  // Determine the energy of the impact.
  background->nextImpact.energy = photon_energy(&background->spectrum, rmf);

  // Determine the time of the background event (impact).
  background->nextImpact.time += rndexp(1./(background->rate));;
  
  return(EXIT_SUCCESS);
}

