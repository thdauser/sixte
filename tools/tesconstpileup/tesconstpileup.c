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


   Copyright 2014 Philippe Peille, IRAP
*/



#include "tesconstpileup.h"

//////////////////////////////////
/** main procedure */
int tesconstpileup_main(){
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  
  // Detector data structure.
  AdvDet *det=NULL;
  
  // Pixel impact list
  PixImpFile* plf=NULL;
  
  // Error status.
  int status=EXIT_SUCCESS;
  
  long ii;
  
  // Register HEATOOL:
  set_toolname("tesconstpileup");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at 
       // most be run once).
       
    // Read parameters using PIL library.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");
    
    // Load detector information
    det=loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    
    // Determine the event list output file.
    char piximplist_filename[MAXFILENAME];
    strcpy(piximplist_filename, par.PixImpList);
    
    // Create origin impactlist_filename
    char impactlist_filename[]="tesconstpileup";
    
    // Open the output file
    plf=openNewPixImpFile(piximplist_filename,
			  par.telescop,
			  par.instrume,
			  par.filter,
			  par.ancrfile,
			  par.respfile,
			  par.XMLFile,
			  impactlist_filename,
			  par.mjdref,
			  par.timezero,
			  0,
			  par.tstop,
			  par.clobber,
			  &status);
    CHECK_STATUS_BREAK(status);
    
    // necessary initialization to have the correct first impact (10 samples after 0)
    double time=(par.preBufferSize-par.triggerSize-90+par.pulseDistance)/det->SampleFreq+.5/det->SampleFreq; 
    ii=0;
    PixImpact piximp;
    while(time<=par.tstop){
      //////////////////////////////////////
      //First impact
      //////////////////////////////////////
      time=time+(par.triggerSize+100-par.pulseDistance)/det->SampleFreq;
      // populate impact structure
      piximp.pixID=0; //fixed pixID (no need to randomize here)
      piximp.time=time;
      piximp.energy=(float)par.energy;
      piximp.ph_id=2*ii;
      piximp.src_id=0;
      piximp.pixposition.x=(-0.5)*det->pix[0].width;
      piximp.pixposition.y=(-0.5)*det->pix[0].height;
      piximp.detposition.x=piximp.pixposition.x+det->pix[0].sx;
      piximp.detposition.y=piximp.pixposition.y+det->pix[0].sy;
      // save impact      
      addImpact2PixImpFile(plf, &piximp, &status);

      //////////////////////////////////////
      //Pileup impact
      //////////////////////////////////////
      time=time+par.pulseDistance/det->SampleFreq;
      piximp.pixID=0;
      piximp.time=time;
      piximp.energy=(float)par.energy;
      piximp.ph_id=2*ii+1;
      piximp.src_id=0;
      piximp.pixposition.x=(-0.5)*det->pix[0].width;
      piximp.pixposition.y=(-0.5)*det->pix[0].height;
      piximp.detposition.x=piximp.pixposition.x+det->pix[0].sx;
      piximp.detposition.y=piximp.pixposition.y+det->pix[0].sy;
      // save impact      
      addImpact2PixImpFile(plf, &piximp, &status);

      ii++;
    }
    
  } while(0); // END of the error handling loop.
  
  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");
  freePixImpFile(&plf, &status);
 
  destroyAdvDet(&det);
  
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

  status=ape_trad_query_string("PixImpList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pixel impact list");
    return(status);
  } 
  strcpy(par->PixImpList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("Telescop", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the telescope");
    return(status);
  } 
  strcpy(par->telescop, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("Instrume", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the instrument");
    return(status);
  } 
  strcpy(par->instrume, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("Filter", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the filter");
    return(status);
  } 
  strcpy(par->filter, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("ancrfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the ancrfile");
    return(status);
  } 
  strcpy(par->ancrfile, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("respfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the respfile");
    return(status);
  } 
  strcpy(par->respfile, sbuffer);
  free(sbuffer); 
  
  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }
  
  status=ape_trad_query_double("tstop", &par->tstop);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstop parameter");
    return(status);
  }
  
  status=ape_trad_query_double("timezero", &par->timezero);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the timezero parameter");
    return(status);
  }
  
  status=ape_trad_query_double("mjdref", &par->mjdref);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the mjdref parameter");
    return(status);
  }
  
  status=ape_trad_query_double("energy", &par->energy);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the energy parameter");
    return(status);
  }
  
  status=ape_trad_query_int("pulseDistance", &par->pulseDistance);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pulse distance");
    return(status);
  }
  
  status=ape_trad_query_int("PreBufferSize", &par->preBufferSize);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pre-buffer size");
    return(status);
  }
  
  status=ape_trad_query_int("TriggerSize", &par->triggerSize);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the trigger size");
    return(status);
  }
  
  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }
  
  return(status);
}
