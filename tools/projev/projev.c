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
*/

#include "sixt.h"
#include "attitude.h"
#include "eventfile.h"
#include "geninst.h"
#include "phproj.h"

#define TOOLSUB projev_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  double MJDREF;
  double TSTART;
  double Exposure;

  int Seed;
  
  char clobber;
};


int projev_getpar(struct Parameters *par);


////////////////////////////////////
/** Main procedure. */
int projev_main() {
  struct Parameters par; // Program parameters

  // Event file.
  EventFile* elf=NULL;

  // Attitude catalog.
  Attitude* ac=NULL;

  // Instrument data structure (containing the focal length, FoV, ...).
  GenInst* inst=NULL;

  // Error status.
  int status=EXIT_SUCCESS;   


  // Register HEATOOL:
  set_toolname("projev");
  set_toolversion("0.04");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=projev_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Initialize the random number generator.
    unsigned int seed=getSeed(par.Seed);
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile, 
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, seed, &status);
    CHECK_STATUS_BREAK(status);

    // Set up the Attitude.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a simple pointing attitude.
      ac=getPointingAttitude(par.MJDREF, par.TSTART, par.TSTART+par.Exposure,
			     par.RA*M_PI/180., par.Dec*M_PI/180., &status);
      CHECK_STATUS_BREAK(status);

    } else {
      // Load the attitude from the given file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      checkAttitudeTimeCoverage(ac, par.MJDREF, par.TSTART, 
				par.TSTART+par.Exposure, &status);
      CHECK_STATUS_BREAK(status);
    }
    // END of setting up the attitude.

    // Set the event file.
    elf=openEventFile(par.EventList, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---


    // --- Beginning of the Sky Projection Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start sky projection process ...\n");

    // Run the pattern projection.
    phproj(inst, ac, elf, par.TSTART, par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Clean up the random number generator.
  sixt_destroy_rng();

  // Destroy the GenInst data structure.
  destroyGenInst(&inst, &status);
  
  // Close the files.
  freeEventFile(&elf, &status);

  // Release memory of Attitude.
  freeAttitude(&ac);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int projev_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mission", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the mission");
    return(status);
  } 
  strcpy(par->Mission, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Instrument", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the instrument");
    return(status);
  } 
  strcpy(par->Instrument, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mode", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the instrument mode");
    return(status);
  } 
  strcpy(par->Mode, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the attitude");
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the right ascension of the telescope "
		   "pointing");
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the declination of the telescope pointing");
    return(status);
  } 

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading MJDREF");
    return(status);
  } 

  status=ape_trad_query_double("TSTART", &par->TSTART);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading TSTART");
    return(status);
  } 

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the exposure time");
    return(status);
  } 

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the seed for the random number generator");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}



