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

#include "gendetsim.h"


////////////////////////////////////
/** Main procedure. */
int gendetsim_main() {
  const double timezero=0.0;

  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Instrument data structure (containing the pixel array, its width, ...).
  GenInst* inst=NULL;

  // Input impact list.
  ImpactFile* ilf=NULL;

  // Output event list file.
  EventFile* elf=NULL;

  // GTI collection.
  GTI* gti=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("gendetsim");
  set_toolversion("0.05");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Determine the appropriate instrument XML definition file.
    char xml_filename[MAXFILENAME];
    sixt_get_XMLFile(xml_filename, par.XMLFile,
		     par.Mission, par.Instrument, par.Mode,
		     &status);
    CHECK_STATUS_BREAK(status);

    // Determine the impact list file.
    char impactlist_filename[MAXFILENAME];
    strcpy(impactlist_filename, par.ImpactList);
    
    // Determine the event list output file.
    char eventlist_filename[MAXFILENAME];
    strcpy(eventlist_filename, par.EventList);

    // Initialize the random number generator.
    unsigned int seed=getSeed(par.Seed);
    sixt_init_rng(seed, &status);
    CHECK_STATUS_BREAK(status);
    
    // Load the instrument configuration.
    inst=loadGenInst(xml_filename, seed, &status);
    CHECK_STATUS_BREAK(status);

    // Use the background if available.
    setGenDetIgnoreBkg(inst->det, 0);

    // Set the start time for the simulation.
    setGenDetStartTime(inst->det, par.TSTART);
    
    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(3, "start detection process ...\n");

    // Open the FITS file with the input impact list:
    ilf=openImpactFile(impactlist_filename, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output event file.
    char telescop[MAXMSG]={""};
    char instrume[MAXMSG]={""};
    if (NULL!=inst->telescop) {
      strcpy(telescop, inst->telescop);
    }
    if (NULL!=inst->instrume) {
      strcpy(instrume, inst->instrume);
    }
    elf=openNewEventFile(eventlist_filename, 
			 telescop, instrume, "Normal", 
			 inst->tel->arf_filename, inst->det->rmf_filename,
			 par.MJDREF, 0.0, par.TSTART, par.TSTART+par.Exposure,
			 inst->det->pixgrid->xwidth, 
			 inst->det->pixgrid->ywidth, 
			 par.clobber, &status);
    CHECK_STATUS_BREAK(status);

    // Set FITS header keywords.
    char keystr[MAXMSG];
    long value;
    sprintf(keystr, "TLMIN%d", elf->cpi);
    value=inst->det->rmf->FirstChannel;
    fits_update_key(elf->fptr, TLONG, keystr, &value, "", &status);
    sprintf(keystr, "TLMAX%d", elf->cpi);
    value=inst->det->rmf->FirstChannel+inst->det->rmf->NumberChannels-1;
    fits_update_key(elf->fptr, TLONG, keystr, &value, "", &status);
    CHECK_STATUS_BREAK(status);
    
    // Store the number of simulated input photons in the FITS header
    // of the output event file.
    fits_update_key(elf->fptr, TLONG, "NPHOTONS", 
		    &ilf->nrows, "number of input photons", 
		    &status);
    CHECK_STATUS_BREAK(status);

    // Store the event type in the FITS header.
    fits_update_key(elf->fptr, TSTRING, "EVTYPE", "PIXEL",
		    "event type", &status);
    CHECK_STATUS_BREAK(status);

    // Define the event list file as output file.
    setGenDetEventFile(inst->det, elf);

    // Loop over all impacts in the FITS file.
    while (ilf->row<ilf->nrows) {

      Impact impact;
      getNextImpactFromFile(ilf, &impact, &status);
      CHECK_STATUS_BREAK(status);

      // Check whether we are still within the requested time interval.
      if (impact.time < par.TSTART) continue;
      if (impact.time > par.TSTART+par.Exposure) break;

      // Photon detection.
      phdetGenDet(inst->det, &impact, par.TSTART+par.Exposure, &status);
      CHECK_STATUS_BREAK(status);

    };
    CHECK_STATUS_BREAK(status);
    // End of loop over all impacts in the input file.

    // Finalize the photon detection.
    phdetGenDet(inst->det, NULL, par.TSTART+par.Exposure, &status);
    CHECK_STATUS_BREAK(status);

    // Store the GTI in the event file.
    gti=newGTI(&status);
    CHECK_STATUS_BREAK(status);
    gti->mjdref=par.MJDREF;
    gti->timezero=timezero;
    appendGTI(gti, par.TSTART, par.TSTART+par.Exposure, &status);
    CHECK_STATUS_BREAK(status);
    saveGTIExt(elf->fptr, "STDGTI", gti, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Clean up the random number generator.
  sixt_destroy_rng();

  freeEventFile(&elf, &status);
  freeImpactFile(&ilf, &status);

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

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("ImpactList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the impact list");
    return(status);
  } 
  strcpy(par->ImpactList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list");
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


