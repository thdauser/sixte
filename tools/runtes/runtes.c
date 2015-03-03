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


   Copyright 2014 Jelle de Plaa, SRON, Thorsten Brand, FAU, Philippe Peille, IRAP
*/

#include "runtes.h"

////////////////////////////////////
/** Main procedure. */
int runtes_main() {

  // Containing all programm parameters read by PIL.
  struct Parameters partmp;
  TESGeneralParameters par;

  // Error status.
  int status=EXIT_SUCCESS;

  int Nstreams=0;
  TESInitStruct* init=NULL;
  TESFitsStream** fitsstream=NULL;
  TESDataStream* stream=NULL;
  fitsfile *ofptr=NULL;
  ReconstructInit* reconstruct_init = NULL;

  int ismonoc=0;
  float monoen=0.;

  int ii;

  // Register HEATOOL:
  set_toolname("runtes");
  set_toolversion("0.05");

  do { // Beginning of the ERROR handling loop (will at
       // most be run once).

    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&partmp);
    CHECK_STATUS_BREAK(status);

    // Copy parameters in general parameters structure
    copyParams2GeneralStruct(partmp,&par);

    // Initialization
    init = newInitStruct(&status);
    CHECK_STATUS_BREAK(status);
    tesinitialization(init,&par,&status);
    CHECK_STATUS_BREAK(status);

    // Generate the data
    stream=newTESDataStream(&status);
    CHECK_STATUS_BREAK(status);

    getTESDataStream(stream,
		     init->impfile,
		     init->profiles,
		     init->det,
		     init->tstart,
		     init->tstop,
		     init->det->npix,
		     par.Nactive,
		     init->activearray,
		     init->Nevts,
		     &ismonoc,
		     &monoen,
		     par.seed,
		     &status);
    CHECK_STATUS_BREAK(status);

    //Write Stream file if asked
    if (par.writeStreamFile) {
      // Get output data arrays
      Nstreams=(par.Nactive+TESFITSMAXPIX-1)/TESFITSMAXPIX;
      fitsstream=(TESFitsStream**)malloc(Nstreams*sizeof(TESFitsStream*));
      if(fitsstream==NULL){
	status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for FITS streams failed");
	CHECK_STATUS_BREAK(status);
      }
      for(ii=0; ii<Nstreams; ii++){
	fitsstream[ii]=newTESFitsStream(&status);
	CHECK_STATUS_BREAK(status);
      }
      createTESFitsStreamFile(&ofptr,
			      par.streamname,
			      init->telescop,
			      init->instrume,
			      init->filter,
			      init->ancrfile,
			      init->respfile,
			      par.XMLFile,
			      par.PixImpList,
			      init->mjdref,
			      init->timezero,
			      init->tstart,
			      init->tstop,
			      par.clobber,
			      &status);
      CHECK_STATUS_BREAK(status);

      // Write the output file

      int restpix=par.Nactive;
      for(ii=0; ii<Nstreams; ii++){
	// print extension name
	sprintf(fitsstream[ii]->name, "ADC%03d", ii+1);

	// calculate number of pixels in the extension
	int extpix=TESFITSMAXPIX;
	if(extpix>restpix){
	  extpix=restpix;
	}

	// allocate memory for fits stream
	allocateTESFitsStream(fitsstream[ii],
			      stream->Ntime,
			      extpix,
			      &status);
	CHECK_STATUS_BREAK(status);

	// transfer values into fits stream
	int pp, tt, ll, id;
	// time array
	for(tt=0; tt<stream->Ntime; tt++){
	  fitsstream[ii]->time[tt]=stream->time[tt];
	}
	// find the ofiginal pixel id (reverse of the activearray)
	for(pp=0; pp<extpix; pp++){
	  id=-1;
	  for(ll=0; ll<init->det->npix; ll++){
	    if(init->activearray[ll]==par.Nactive-restpix+pp){
	      id=ll;
	    }
	  }
	  if(id==-1){
	    status=EXIT_FAILURE;
	    SIXT_ERROR("Pixel ID not found.");
	    CHECK_STATUS_BREAK(status);
	  }
	  CHECK_STATUS_BREAK(status);
	  fitsstream[ii]->pixID[pp]=id;

	  // adc_value array
	  for(tt=0; tt<stream->Ntime; tt++){
	    fitsstream[ii]->adc_value[pp][tt]=stream->adc_value[tt][par.Nactive-restpix+pp];
	  }
	}
	CHECK_STATUS_BREAK(status);

	double timeres=1./init->det->SampleFreq;
	writeTESFitsStream(ofptr,
			   fitsstream[ii],
			   init->tstart,
			   init->tstop,
			   timeres,
			   init->Nevts,
			   ismonoc,
			   monoen,
			   &status);
	CHECK_STATUS_BREAK(status);
	restpix=restpix-TESFITSMAXPIX;
      }
      CHECK_STATUS_BREAK(status);

      fits_close_file(ofptr, &status);
      CHECK_STATUS_BREAK(status);
    }

    //Initialize reconstruction if necessary
    if (partmp.Reconstruct){
		reconstruct_init = newReconstructInit(&status);
		CHECK_STATUS_BREAK(status);
		initializeReconstruction(reconstruct_init,partmp.OptimalFilterFile,partmp.PulseLength,
				partmp.PulseTemplateFile,partmp.Threshold,partmp.Calfac,partmp.NormalExclusion,
	    		partmp.DerivateExclusion,partmp.SaturationValue,&status);
		CHECK_STATUS_BREAK(status);
    }

    //Trigerring part
    triggerWithImpact(stream,&par,init,monoen,partmp.WriteRecordFile,partmp.Reconstruct,
    		reconstruct_init,partmp.TesEventFile,partmp.EventListSize,partmp.Identify,&status);

  } while(0); // END of the error handling loop.

  freeTESInitStruct(&init,&status);
  freeReconstructInit(reconstruct_init);
  destroyTESDataStream(stream);
  for(ii=0; ii<Nstreams; ii++){
    destroyTESFitsStream(fitsstream[ii]);
  }
  free(fitsstream);

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

  status=ape_trad_query_bool("Reconstruct", &par->Reconstruct);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the Reconstruct parameter");
    return(status);
  }

  status=ape_trad_query_bool("WriteRecordFile", &par->WriteRecordFile);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the WriteRecordFile parameter");
    return(status);
  }

  if (par->WriteRecordFile){
	  status=ape_trad_query_string("TesTriggerfile", &sbuffer);
	  if (EXIT_SUCCESS!=status) {
		  SIXT_ERROR("failed reading the name of the TES Trigger output file");
		  return(status);
	  }
	  strcpy(par->tesTriggerFile, sbuffer);
	  free(sbuffer);
  }

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

  status=ape_trad_query_bool("WriteStreamFile", &par->writeStreamFile);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the WriteStreamFile parameter");
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

  if (par->Reconstruct) {
	  status=ape_trad_query_string("TesEventFile", &sbuffer);
	  if (EXIT_SUCCESS!=status) {
		  SIXT_ERROR("failed reading the name of the event file");
		  return(status);
	  }
	  strcpy(par->TesEventFile, sbuffer);
	  free(sbuffer);

	  status=ape_trad_query_string("OptimalFilterFile", &sbuffer);
	  if (EXIT_SUCCESS!=status) {
		  SIXT_ERROR("failed reading the name of the optimal filter file");
		  return(status);
	  }
	  strcpy(par->OptimalFilterFile, sbuffer);
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

	  status=ape_trad_query_bool("Identify", &par->Identify);
	  if (EXIT_SUCCESS!=status) {
		  SIXT_ERROR("failed reading the Identify parameter");
		  return(status);
	  }
  }

  if (par->WriteRecordFile==0 && par->Reconstruct==0){
	  SIXT_ERROR("You have to at least choose to save the records or to reconstruct the events. Aborting simulation");
	  return(EXIT_FAILURE);
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

  int seed;
  status=ape_trad_query_int("Seed", &seed);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the Seed");
    return(status);
  }
  if(seed==-1){
    // Initialize with system time
    struct timeval tv;
    gettimeofday(&tv,NULL);
    par->seed=1000000*tv.tv_sec+tv.tv_usec;
  }else{
    par->seed=(unsigned long int)seed;
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

  par->writeStreamFile=partmp.writeStreamFile;
  par->clobber=partmp.clobber;
  par->history=partmp.history;

  par->seed=partmp.seed;
}

