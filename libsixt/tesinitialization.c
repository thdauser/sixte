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

#include "tesinitialization.h"

/** Initializes the different variables necessary fo the simulations. Depending
    on the tool calling this function, not all the variables are set. */
void tesinitialization(TESInitStruct* const init,TESGeneralParameters* const par, int* const status){
  int ii;

  // Open the pixel impact file
  init->impfile=openPixImpFile(par->PixImpList, READONLY,status);
  CHECK_STATUS_VOID(*status);
  // Read keywords from input file
  sixt_read_fits_stdkeywords_obsolete(init->impfile->fptr,
			     init->telescop,
			     init->instrume,
			     init->filter,
			     init->ancrfile,
			     init->respfile,
			     &(init->mjdref),
			     &(init->timezero),
			     &(init->tstart),
			     &(init->tstop), 
			     status);
  CHECK_STATUS_VOID(*status);
  if(par->check_times){
	  printf("Pixel impact file reaches from %lfs-%lfs .\n", init->tstart, init->tstop);
	  if(init->tstart>par->tstart){
		  puts("Program parameter tstart smaller than in input file.");
	  }else{
		  init->tstart=par->tstart;
	  }
	  if(init->tstop<par->tstop){
		  puts("Program parameter tstop larger than in input file.");
	  }else{
		  init->tstop=par->tstop;
	  }
	  printf("Simulate from %lfs-%lfs .\n", init->tstart, init->tstop);
  } else {
	  init->tstart = par->tstart;
	  init->tstop = par->tstop;
  }
  
  // Load the detector structure
  init->det=loadAdvDet(par->XMLFile,status);
  CHECK_STATUS_VOID(*status);
    
  // construct array of active pixels
  init->activearray=(int*)malloc(init->det->npix*sizeof(int));
  if(init->activearray==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for array of actove pixels failed");
    CHECK_STATUS_VOID(status);
  }
  if(par->Nactive==-1){
    par->Nactive=init->det->npix;
    par->nlo=0;
    par->nhi=init->det->npix-1;
  }
  int act=0;
  for(ii=0; ii<init->det->npix; ii++){
    if(ii<par->nlo || ii>par->nhi){
      init->activearray[ii]=-1;
    } else {
      init->activearray[ii]=act;
      act++;
    }
  }
  if(act!=par->Nactive){
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Number of active pixels has been corrupted. Counted %d instead of %d pixels.", act, par->Nactive);
    SIXT_ERROR(msg);
    CHECK_STATUS_VOID(*status);
  }
    
  // construct array for event numbers
  init->Nevts=(long*)malloc(init->det->npix*sizeof(long));
  if(init->Nevts==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for array of event numbers failed");
    CHECK_STATUS_VOID(*status);
  }
  for(ii=0; ii<init->det->npix; ii++){
    init->Nevts[ii]=0;
  }
    
  init->profiles=newTESProfiles(status);
  CHECK_STATUS_VOID(*status);
    
  for(ii=0; ii<init->det->npix; ii++){
    // Test if profile is already loaded, if not, load it
    int versionindex=findTESProfileVersionIndex(init->profiles, init->det->pix[ii].version);
    if(init->activearray[ii]>-1 && versionindex<0){
      char profilename[MAXFILENAME];
      sprintf(profilename, "%s%s", init->det->filepath, init->det->tesproffilename);
      readTESProfiles(profilename,
		      init->det->pix[ii].version, 
		      init->profiles, 
		      status);
      CHECK_STATUS_VOID(*status);
      versionindex=findTESProfileVersionIndex(init->profiles, init->det->pix[ii].version);
      if(versionindex<0){
	*status=EXIT_FAILURE;
	SIXT_ERROR("New profile not loaded.");
	CHECK_STATUS_VOID(*status);
      }
      init->det->pix[ii].profVersionID=versionindex;
    } else {
      init->det->pix[ii].profVersionID=versionindex;
    }
  }
  if(*status!=EXIT_SUCCESS){
    SIXT_ERROR("Failed reading pulse profile templates.");
    CHECK_STATUS_VOID(*status);
  }
    
  //Open output files if needed
  SixtStdKeywords* keywords = buildSixtStdKeywords(init->telescop,init->instrume,init->filter,
		  init->ancrfile,init->respfile,"NONE",init->mjdref,init->timezero,init->tstart,init->tstop,status);
  if(par->WriteRecordFile){
	  init->record_file= opennewTesTriggerFile(par->tesTriggerFile,keywords,par->XMLFile,par->PixImpList,
			  par->triggerSize,par->preBufferSize,init->det->SampleFreq,par->clobber,status);
	  CHECK_STATUS_VOID(*status);
  }
  if (par->Reconstruct){
	  init->event_file = opennewTesEventFile(par->TesEventFile,keywords,par->clobber,status);
	  CHECK_STATUS_VOID(*status);
  }
}

/** Constructor. Returns a pointer to an empty TESInitStruct data
    structure. */
TESInitStruct* newInitStruct(int* const status){
  TESInitStruct* init = (TESInitStruct*)malloc(sizeof(TESInitStruct));
  if (NULL==init) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TESInitStruct failed");
    return(init);
  }

  // Initialize pointers with NULL.
  init->impfile       =NULL;
  init->profiles      =NULL;
  init->det           =NULL;
  init->activearray   =NULL;
  init->Nevts         =NULL;
  init->record_file   =NULL;
  init->event_file    =NULL;
  
  // Initialize values.
  init->mjdref	   =0;
  init->timezero   =0;
  init->tstart	   =0;
  init->tstop      =0;

  return(init);
   
}

/** Destructor. */
void freeTESInitStruct(TESInitStruct** const init, int* const status){  
  if (NULL!=*init) {
    //Free pointers
    freePixImpFile(&((*init)->impfile), status);
    destroyTESProfiles((*init)->profiles);
    destroyAdvDet(&((*init)->det));
    freeTesTriggerFile(&((*init)->record_file), status);
    freeTesEventFile((*init)->event_file,status);
    
    if((*init)->activearray!=NULL){
      free((*init)->activearray);
      (*init)->activearray=NULL;
    }
    
    if((*init)->Nevts!=NULL){
      free((*init)->Nevts);
      (*init)->Nevts=NULL;
    }

    // Erase values.
    (*init)->mjdref     =0;
    (*init)->timezero   =0;
    (*init)->tstart     =0;
    (*init)->tstop      =0;
  }
  free(*init);
  *init=NULL;
}
