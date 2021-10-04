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

#include "sixt.h"
#include "attitude.h"
#include "eventfile.h"
#include "geninst.h"
#include "phproj.h"

#define TOOLSUB projev_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char RawData[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char AdvXml[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec, rollangle;

  double MJDREF;
  double TSTART;
  double Exposure;

  int Seed;

  char clobber;
  char ProjCenter;
};


int projev_getpar(struct Parameters *par);


////////////////////////////////////
/** Main procedure. */
int projev_main() {
  struct Parameters par; // Program parameters

  // Event file.
  EventFile* elf=NULL;

  // TES event file.
  TesEventFile* tes_elf=NULL;

  // Attitude catalog.
  Attitude* ac=NULL;

  // Instrument data structure (containing the focal length, FoV, ...).
  GenInst* inst=NULL;

  // Advanced detectector data structure
  AdvDet* adv_det=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("projev");
  set_toolversion("0.05");


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

    // Determine the instrument XML definition file.
    char xml_filename[MAXFILENAME];
    strcpy(xml_filename, par.XMLFile);

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
			     par.RA*M_PI/180., par.Dec*M_PI/180., par.rollangle*M_PI/180., &status);
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
    strcpy(ucase_buffer,par.AdvXml);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")){
    	// Load classic event file
    	elf=openEventFile(par.RawData, READWRITE, &status);
    	if(status==COL_NOT_FOUND){
    		SIXT_WARNING("You may have given this tool a TES event file as input without giving it an advanced xml file");
    	}
    } else {
    	// Load Tes event file
    	tes_elf=openTesEventFile(par.RawData,READWRITE,&status);
    }
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---


    // --- Beginning of the Sky Projection Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start sky projection process ...\n");

    // Run the pattern projection.
    if (NULL!=elf){
    	phproj(inst, ac, elf, par.TSTART, par.Exposure,&status);
    } else {
    	adv_det = loadAdvDet(par.AdvXml,&status);
		CHECK_STATUS_BREAK(status);
    	phproj_advdet(inst,adv_det,ac,tes_elf,par.TSTART,par.Exposure,par.ProjCenter,&status);
    }
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Clean up the random number generator.
  sixt_destroy_rng();

  // Destroy the GenInst data structure.
  destroyGenInst(&inst, &status);

  // Destroy the AdvDet data structure
  destroyAdvDet(&adv_det);

  // Close the files.
  freeEventFile(&elf, &status);
  freeTesEventFile(tes_elf,&status);

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

  // check if any obsolete keywords are given
  sixt_check_obsolete_keyword(&status);
  CHECK_STATUS_RET(status,EXIT_FAILURE);


  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("RawData", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  }
  strcpy(par->RawData, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  }
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("AdvXml", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the advanced detector XML file");
    return(status);
  }
  strcpy(par->AdvXml, sbuffer);
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

  query_simput_parameter_float("rollangle",&(par->rollangle),&status);


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

  status=ape_trad_query_bool("ProjCenter", &par->ProjCenter);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the ProjCenter parameter");
	  return(status);
  }


  return(status);

}
