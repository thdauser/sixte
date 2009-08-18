#include "pattern_recombination.h"


int pattern_recombination_main() {
  struct Parameters parameters;

  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("pattern_recombination");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = pattern_recombination_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    
  } while(0); // End of error handling loop


  // --- Clean Up ---

  return(status);
}



int pattern_recombination_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input file!\n", status);
  }

  return(status);
}





