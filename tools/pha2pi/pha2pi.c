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

#include "pha2pi.h"


int pha2pi_main()
{
  // Containing all program parameters read by PIL
  struct Parameters par;

  // Pha2PI file
  Pha2Pi* p2p=NULL;

  // Input event file.
  EventFile* evtfile=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("pha2pi");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    pha2pi_getpar(&par,&status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    unsigned int seed=getSeed(par.Seed);

    // Determine the Pha2Pi file.
    char pha2pi_filename[MAXFILENAME];
    strcpy(pha2pi_filename, par.Pha2Pi);

    // Determine the output file.
    char evtfile_filename[MAXFILENAME];
    strcpy(evtfile_filename, par.EvtFile);

    // Determine the Pha2Pi file.
    char RSPPath[MAXFILENAME];
    strcpy(RSPPath, par.RSPPath);

    // Determine the Pha2Pi file.
    char RESPfile[MAXFILENAME];
    strcpy(RESPfile, par.RESPfile);

    headas_chat(3, "start PHA to PI correction ...\n");

    // LOAD Pha2PI Correction File if existent
    p2p = initPha2Pi(pha2pi_filename, seed, &status);
    CHECK_STATUS_BREAK_WITH_FITSERROR(status);

    // Open the input event file for read and write.
    evtfile=openEventFile(evtfile_filename, READWRITE, &status);
    CHECK_STATUS_BREAK_WITH_FITSERROR(status);

    // Run PI correction on evtfile.
    pha2pi_correct_eventfile( evtfile, p2p, RSPPath, RESPfile, &status);
    CHECK_STATUS_BREAK(status);


  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&evtfile, &status);
  freePha2Pi(&p2p);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


void pha2pi_getpar(struct Parameters* const par, int* const status){
  query_simput_parameter_file_name_buffer("EvtFile", par->EvtFile, MAXFILENAME, status);
  query_simput_parameter_file_name_buffer("Pha2Pi", par->Pha2Pi, MAXFILENAME, status);
  query_simput_parameter_file_name_buffer("RSPPath", par->RSPPath, MAXFILENAME, status);
  query_simput_parameter_file_name_buffer("RESPfile", par->RESPfile, MAXFILENAME, status);
  query_simput_parameter_int("Seed", &par->Seed, status );
  query_simput_parameter_bool("clobber", &par->clobber, status);
}
