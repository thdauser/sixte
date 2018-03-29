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

#include "evpat.h"


int evpat_main() 
{
  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Input event file.
  EventFile* elf=NULL;

  // Output event pattern file.
  EventFile* plf=NULL;

  // Instrument.
  GenInst* inst=NULL;

  // GTI collection.
  GTI* gti=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("evpat");
  set_toolversion("0.06");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile, 
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Load the instrument configuration.
    unsigned int seed=getSeed(par.Seed);
    inst=loadGenInst(xml_filename, seed, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the event list file name.
    char rawdata_filename[MAXFILENAME];
    strcpy(rawdata_filename, par.RawData);

    // Determine the output file.
    char evtfile_filename[MAXFILENAME];
    strcpy(evtfile_filename, par.EvtFile);

    // Determine the output file.
    char pha2pi_filename[MAXFILENAME];
    strcpy(pha2pi_filename, par.Pha2Pi);


    // Load the GTI extension from the event file.
    gti=loadGTI(rawdata_filename, &status);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "start pattern recombination ...\n");

    // Open the input event file.
    elf=openEventFile(rawdata_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read the timing keywords.
    char comment[MAXMSG];
    double mjdref, timezero, tstart, tstop;
    fits_read_key(elf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output file.
    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    plf=openNewEventFile(evtfile_filename,
			 telescop, instrume, "Normal",
			 inst->tel->arf_filename,
			 inst->det->rmf_filename,
			 mjdref, timezero, tstart, tstop,			   
			 inst->det->pixgrid->xwidth,
			 inst->det->pixgrid->ywidth,
			 par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Pattern recombination.
    phpat(inst->det, elf, plf, pha2pi_filename, par.SkipInvalids, &status);
    CHECK_STATUS_BREAK(status);

    // Store the GTI in the pattern file.
    saveGTIExt(plf->fptr, "STDGTI", gti, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&elf, &status);
  freeEventFile(&plf, &status);
 
  destroyGenInst(&inst, &status);
  freeGTI(&gti);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // check if any obsolete keywords are given
  sixt_check_obsolete_keyword(&status);
  CHECK_STATUS_RET(status,EXIT_FAILURE);

//  status=ape_trad_query_file_name("RawData", &sbuffer);
//  if (EXIT_SUCCESS!=status) {
//    SIXT_ERROR("failed reading the name of the input event list");
//    return(status);
//  }
//  strcpy(par->RawData, sbuffer);
//  free(sbuffer);
//
//  status=ape_trad_query_file_name("EvtFile", &sbuffer);
//  if (EXIT_SUCCESS!=status) {
//    SIXT_ERROR("failed reading the name of the output pattern list");
//    return(status);
//  }
//  strcpy(par->EvtFile, sbuffer);
//  free(sbuffer);

  query_simput_parameter_file_name_buffer("RawData", par->RawData, MAXFILENAME, &status);
  query_simput_parameter_file_name_buffer("EvtFile", par->EvtFile, MAXFILENAME, &status);
  query_simput_parameter_file_name_buffer("Pha2Pi", par->Pha2Pi, MAXFILENAME, &status);

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

  status=ape_trad_query_bool("SkipInvalids", &par->SkipInvalids);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the SkipInvalids parameter");
    return(status);
  }

  status=ape_trad_query_int("Seed", &par->Seed);
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


