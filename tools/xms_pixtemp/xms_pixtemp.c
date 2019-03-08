/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "xms_pixtemp.h"


int xms_pixtemp_main() {
  struct Parameters parameters;
  EventFile* elf=NULL;
  FILE* output_file=NULL;

  int status=EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("xms_pixtemp");
  set_toolversion("0.02");


  do { // ERROR handling loop

    // Read parameters by PIL:
    status=xms_pixtemp_getpar(&parameters);
    CHECK_STATUS_BREAK(status);

    // Initialize the random number generator.
    sixt_init_rng((int)time(NULL), &status);
    CHECK_STATUS_BREAK(status);

    // Open the event file.
    elf=openEventFile(parameters.EventList, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // Read the EBOUNDS from the detector response file.
    struct RMF* rmf=loadRMF(parameters.RSP, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output file.
    output_file=fopen(parameters.OutputFile, "w+");
    if (NULL==output_file) {
      status=EXIT_FAILURE;
      SIXT_ERROR("opening the output file failed");
      break;
    }

    // Loop over all events in the event file.
    long row;
    for (row=0; row<elf->nrows; row++) {

      // Read the next event from the FITS file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      /*      if ((1==event.array) && // Only events from the inner array.
	      (event.xi == parameters.pixx) && (event.xi == parameters.pixx)) { */
      fprintf(output_file, " %lf\t%lf\n", event.time,
	      getEBOUNDSEnergy(event.pha, rmf, &status));
      CHECK_STATUS_BREAK(status);
      /* } */

    } // End of loop over all events in the event file

  } while(0); // End of error handling loop


  // --- Clean Up ---

  // Close the event file.
  freeEventFile(&elf, &status);

  // Close the output file.
  if (NULL!=output_file) {
    fclose(output_file);
    output_file=NULL;
  }

  // Clean up the random number generator.
  sixt_destroy_rng();

  return(status);
}


int xms_pixtemp_getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS;

  if ((status=PILGetFname("EventList", parameters->EventList))) {
    SIXT_ERROR("failed reading the name of the input file");
  }

  else if ((status=PILGetFname("OutputFile", parameters->OutputFile))) {
    SIXT_ERROR("failed reading the name of the output file");
  }

  else if ((status=PILGetFname("RSP", parameters->RSP))) {
    SIXT_ERROR("failed reading the name of the detector response file");
  }

  else if ((status=PILGetInt("pixx", &parameters->pixx))) {
    SIXT_ERROR("failed reading x-coordinate of the pixel to be analysed");
  }

  else if ((status=PILGetInt("pixy", &parameters->pixy))) {
    SIXT_ERROR("failed reading y-coordinate of the pixel to be analysed");
  }

  return(status);
}
