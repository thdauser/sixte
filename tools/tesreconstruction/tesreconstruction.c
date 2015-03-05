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


   Copyright 2015 Philippe Peille, IRAP
*/

#include "tesreconstruction.h"

////////////////////////////////////
/** Main procedure. */
int tesreconstruction_main() {
  // Containing all programm parameters read by PIL.
  struct Parameters par;


  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("tesreconstruction");
  set_toolversion("0.05");

  do { // Beginning of the ERROR handling loop (will at
       // most be run once).
    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Sixt standard keywords structure
    SixtStdKeywords* keywords = newSixtStdKeywords(&status);
    CHECK_STATUS_BREAK(status);

    // Open record file
    TesTriggerFile* record_file = openexistingTesTriggerFile(par.RecordFile,keywords,&status);
    CHECK_STATUS_BREAK(status);

    //Open outfile
    TesEventFile * outfile = opennewTesEventFile(par.TesEventFile,
    		keywords,
    		par.clobber,
    		&status);
    CHECK_STATUS_BREAK(status);

    // Initialize data structures needed for pulse filtering
    ReconstructInit* reconstruct_init = newReconstructInit(&status);
    CHECK_STATUS_BREAK(status);
    initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.PulseLength,
    		par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
    		par.DerivateExclusion,par.SaturationValue,&status);
    CHECK_STATUS_BREAK(status);

    // Build up TesRecord to read the file
    TesRecord* record = newTesRecord(&status);
    allocateTesRecord(record,record_file->trigger_size,record_file->delta_t,0,&status);
    CHECK_STATUS_BREAK(status);

    //Build up TesEventList to recover the results of the reconstruction
    TesEventList* event_list = newTesEventList(&status);
    allocateTesEventListTrigger(event_list,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);


    // Iterate of records and do the reconstruction
    while(getNextRecord(record_file,record,&status)){
    	reconstructRecord(record,event_list,reconstruct_init,0,&status);
    	saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
    	//Reinitialize event list
    	event_list->index=0;
    }
    CHECK_STATUS_BREAK(status);

    //Free memory
    freeReconstructInit(reconstruct_init);
    freeTesTriggerFile(&record_file,&status);
    freeTesEventFile(outfile,&status);
    freeTesEventList(event_list);
    freeTesRecord(&record);
    freeSixtStdKeywords(keywords);
    CHECK_STATUS_BREAK(status);


  } while(0); // END of the error handling loop.


  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(status);
  }
}

int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("OptimalFilterFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the optimal filter file");
    return(status);
  }
  strcpy(par->OptimalFilterFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("RecordFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the optimal filter file");
    return(status);
  }
  strcpy(par->RecordFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("TesEventFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  }
  strcpy(par->TesEventFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("PulseTemplateFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pulse template file");
    return(status);
  }
  strcpy(par->PulseTemplateFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("PulseLength", &par->PulseLength);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the PulseLength parameter");
	  return(status);
  }

  status=ape_trad_query_double("Threshold", &par->Threshold);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the Threshold parameter");
	  return(status);
  }

  status=ape_trad_query_double("Calfac", &par->Calfac);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the Calfac parameter");
	  return(status);
  }

  status=ape_trad_query_int("EventListSize", &par->EventListSize);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the EventListSize parameter");
	  return(status);
  }

  status=ape_trad_query_int("NormalExclusion", &par->NormalExclusion);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the NormalExclusion parameter");
	  return(status);
  }

  status=ape_trad_query_int("DerivateExclusion", &par->DerivateExclusion);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the DerivateExclusion parameter");
	  return(status);
  }

  status=ape_trad_query_double("SaturationValue", &par->SaturationValue);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the SaturationValue parameter");
	  return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the clobber parameter");
	  return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }

  return(status);
}


