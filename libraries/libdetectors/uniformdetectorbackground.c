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

  return(status);
}



int cleanupUniformDetectorBackground(UniformDetectorBackground* background)
{
  // Release allocated memory.
  cleanupSpectrum(&background->spectrum);

  return(EXIT_SUCCESS);
}

