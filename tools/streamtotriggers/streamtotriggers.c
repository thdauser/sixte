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

#include "streamtotriggers.h"
#include "tesinitialization.h"

////////////////////////////////////
/** Main procedure. */
int streamtotriggers_main() {
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  TESGeneralParameters generic_par;
  
  // Error status.
  int status=EXIT_SUCCESS;
  
  //Pointers
  TESDataStream* stream=NULL;
  AdvDet *det=NULL;
  TesStreamFile* tesfile=NULL;
  TESInitStruct* init=NULL;

  float monoen;
  
  int ii;
  int hdu_type=0;
  int act=0;

  char comment[MAXMSG];

  // Register HEATOOL:
  set_toolname("streamtotriggers");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at 
       // most be run once).
    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Load the detector structure
    det=loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    
    // Check and eventually modifies the active pixels
    if(par.Nactive==-1){
      par.Nactive=det->npix;
      par.nlo=0;
      par.nhi=det->npix-1;
    }
    for(ii=0; ii<det->npix; ii++){
      if(ii>=par.nlo && ii<=par.nhi){
	act++;
      }
    }
    if(act!=par.Nactive){
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "Number of active pixels has been corrupted. Counted %d instead of %d pixels.", act, par.Nactive);
      SIXT_ERROR(msg);
      CHECK_STATUS_BREAK(status);
    }
    
    // Initialization structure necessary for trigger step
    init = newInitStruct(&status);
    CHECK_STATUS_BREAK(status);

    //Open TesStream File
    tesfile = openTesStreamFile(par.streamname,par.nlo,par.nhi,READONLY, &status);
    CHECK_STATUS_BREAK(status);


    // Read keywords from tes file
    if (fits_movabs_hdu(tesfile->fptr,1, &hdu_type, &status)) break;
    sixt_read_fits_stdkeywords(tesfile->fptr,
			       init->telescop,
			       init->instrume,
			       init->filter,
			       init->ancrfile,
			       init->respfile,
			       &init->mjdref,
			       &init->timezero,
			       &init->tstart,
			       &init->tstop,
			       &status);
    CHECK_STATUS_BREAK(status);

    //Read monochromatic energy
    if (fits_movabs_hdu(tesfile->fptr,2, &hdu_type, &status)) break;
    fits_read_key(tesfile->fptr, TFLOAT, "MONOEN", &monoen, comment, &status);
    if (status==KEY_NO_EXIST) {
      puts("No monochromatic energy keyword in stream file, defaulting to 0.");
      monoen=0;
      status=EXIT_SUCCESS;
    }
    CHECK_STATUS_BREAK(status);



    double tstart_stream;
    printf("TES stream file reaches from %lfs-%lfs .\n", init->tstart, init->tstop);
    tstart_stream = init->tstart;
    if(init->tstart>par.tstart){
      if(init->tstart>par.tstop){
	SIXT_ERROR("Program parameter tstop smaller than tstart in TES stream file -> abort");
	status=EXIT_FAILURE;
	CHECK_STATUS_BREAK(status);
      }
      puts("Program parameter tstart smaller as in TES stream file.");    
    }else{
      init->tstart=par.tstart;
    }
    if(init->tstop<par.tstop){
      if(init->tstop<par.tstart){
	SIXT_ERROR("Program parameter tstart larger than tstop in TES stream file -> abort");
	status=EXIT_FAILURE;
	CHECK_STATUS_BREAK(status);
      }
      puts("Program parameter tstop larger as in TES stream file.");
    }else{
      init->tstop=par.tstop;
    }
    printf("Retrieving TES Stream from %lfs-%lfs .\n", init->tstart, init->tstop);

    //Create TESDataStream from file
    stream = newTESDataStream(&status);
    CHECK_STATUS_BREAK(status);
    stream = generateTESDataStreamFromFile(stream,tesfile,tstart_stream,init->tstart,init->tstop,
					   det->SampleFreq,&status);
    CHECK_STATUS_BREAK(status);

    // Copy parameters in general parameters structure
    copyParams2GeneralStruct(par,&generic_par);
    triggerWithImpact(stream,&generic_par,init,monoen,1,0,NULL,"dummy",0,&status);

  } while(0); // END of the error handling loop.

  //Free memory
  destroyAdvDet(&det);
  freeTESInitStruct(&init,&status);
  destroyTESDataStream(stream);
  freeTesStreamFile(&tesfile,&status);

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

  status=ape_trad_query_string("Streamfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the Stream file");
    return(status);
  } 
  strcpy(par->streamname, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("TesTriggerfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the TES Trigger output file");
    return(status);
  } 
  strcpy(par->tesTriggerFile, sbuffer);
  free(sbuffer);
  
  status=ape_trad_query_int("TriggerSize", &par->triggerSize);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the TriggerSize parameter");
    return(status);
  }

  status=ape_trad_query_int("PreBufferSize", &par->preBufferSize);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the PreBufferSize parameter");
    return(status);
  }
  
  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  status=ape_trad_query_double("tstart", &par->tstart);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstart parameter");
    return(status);
  }
  
  status=ape_trad_query_double("tstop", &par->tstop);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstop parameter");
    return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }
  
  char *pix=NULL;
  status=ape_trad_query_string("pixels", &pix);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pixel parameter");
    return(status);
  }
  if(!strcmp(pix, "all")){
    par->Nactive=-1;
  } else {
    char *ptr=strchr(pix, '-');
    if(ptr!=NULL){
      sscanf(pix, "%d - %d", &(par->nlo), &(par->nhi));
      par->nlo--;
      par->nhi--;
      par->Nactive=par->nhi-par->nlo+1;
      printf("Activate pixels %d-%d (%d pixels).\n", par->nlo+1, par->nhi+1, par->Nactive);
    } else {
      sscanf(pix, "%d", &(par->nlo));
      par->nlo--;
      par->nhi=par->nlo;
      par->Nactive=1;
    }
    if(par->Nactive<0 || par->nlo<0){
      status=EXIT_FAILURE;
      SIXT_ERROR("Input values for active pixels corrupted.");
      return(status);
    }
  }
  free(pix);

  return(status);
}

/** Copies the parameters contained in the local parameter structure into the
    more general one*/
void copyParams2GeneralStruct(const struct Parameters partmp, TESGeneralParameters* const par){
  strcpy(par->PixImpList,partmp.PixImpList);
  strcpy(par->XMLFile,partmp.XMLFile);
  strcpy(par->streamname,partmp.streamname);
  strcpy(par->tesTriggerFile,partmp.tesTriggerFile);

  strcpy(par->activePixels,partmp.activePixels);
  par->Nactive=partmp.Nactive;
  par->Npix=partmp.Npix;
  par->nlo=partmp.nlo;
  par->nhi=partmp.nhi;
  par->triggerSize=partmp.triggerSize;
  par->preBufferSize=partmp.preBufferSize;

  par->tstart=partmp.tstart;
  par->tstop=partmp.tstop;

  par->clobber=partmp.clobber;
  par->history=partmp.history;
}

